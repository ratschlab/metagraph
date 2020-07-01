#!/usr/bin/env python3

import argparse
import os.path
import subprocess
import time

import googleapiclient.discovery


def count_pending_operations(compute, project, zone, operation_ids):
    id_filter = ' OR '.join(f'(id={id})' for id in operation_ids)
    result = compute.zoneOperations().list(
        project=project,
        zone=zone,
        filter='(status!=DONE) AND ' + id_filter).execute()

    items = result.get('items')
    return items.count if items else 0


def wait_for_all_operations(compute, project, zone, ids):
    if not ids:
        print('There are no operations to wait for to finish. This is probably an error.')
        return
    print('Waiting for *all* operations to finish...', end='', flush=True)
    sleep_time_sec = 1
    while True:
        if count_pending_operations(compute, project, zone, ids) == 0:
            print("done.")
            return
        print('.', end='', flush=True)
        time.sleep(sleep_time_sec)
        if sleep_time_sec < 16:
            sleep_time_sec *= 2


def wait_for_operation(compute, project, zone, operation):
    print('Waiting for operation to finish...', end='', flush=True)
    while True:
        result = compute.zoneOperations().get(
            project=project,
            zone=zone,
            operation=operation
        ).execute()
        print('.', end='', flush=True)
        if result['status'] == 'DONE':
            print("done.")
            if 'error' in result:
                raise Exception(result['error'])
            return result

        time.sleep(1)


def list_instances(compute, project, zone, name):
    """ Lists the currently available instance """
    name_filter = ''
    if name != '':
        name_filter = 'name=' + name
    result = compute.instances().list(project=project, zone=zone, filter=name_filter).execute()
    instances = result['items'] if 'items' in result else None
    if not instances:
        print('No instances found.')
        return
    print(f'Instances in project {project} and zone {zone}:');
    for instance in instances:
        print(' - ' + instance['name'])


def send_delete_request(compute, project, zone, name):
    """ Sends a request to delete the instance with the given name """
    print(f'Deleting instance {name} ...')
    return compute.instances().delete(
        project=project,
        zone=zone,
        instance=name).execute()


def delete_instances(compute, project, zone, name, count):
    ids = []
    for i in range(count):
        operation = send_delete_request(compute, project, zone, name + '-' + str(i))
        ids.append(operation['id'])

    wait_for_all_operations(compute, project, zone, ids)


def send_start_request(compute, project, zone, name):
    """ Sends a request to start the instance with the given name """
    print(f'Starting instance {name} ...')
    return compute.instances().start(
        project=project,
        zone=zone,
        instance=name).execute()


def start_instances(compute, project, zone, name, count):
    """ Starts count instances with the given name prefix. The names will be name-0, ..., name-count-1 """
    ids = []
    for i in range(count):
        operation = send_start_request(compute, project, zone, name + '-' + str(i))
        ids.append(operation['id'])
    wait_for_all_operations(compute, project, zone, ids)


def send_stop_request(compute, project, zone, name):
    """ Sends a request to stop  the instance with the given name """
    print(f'Stopping instance {name} ...')
    return compute.instances().stop(
        project=project,
        zone=zone,
        instance=name).execute()


def stop_instances(compute, project, zone, name, count):
    """ Stops count instances with the given name prefix. The names will be name-0, ..., name-count-1 """
    ids = []
    for i in range(count):
        operation = send_stop_request(compute, project, zone, name + '-' + str(i))
        ids.append(operation['id'])
    wait_for_all_operations(compute, project, zone, ids)


def send_create_request(compute, project, zone, name, script_dir, server_host, num_instances):
    """ Send a request to creates a new instance with the given parameters """
    print(f'Creating instance {name} ...')
    # Get the metagraph Ubuntu 18.04 TLS image
    image_response = compute.images().get(project='metagraph', image='mg-image-40').execute()
    # image_response = compute.snapshots().get(project='metagraph', snapshot='metagraph7').execute()
    source_disk_image = image_response['selfLink']

    # Configure the machine: 1vCPU, 3.75GB of RAM
    machine_type = f"zones/{zone}/machineTypes/n1-standard-1"
    metadata = [
        {
            'key': 'instance_id',
            'value': name.split('-')[-1]
        },
        {
            'key': 'num_instances',
            'value': num_instances
        },
        {
            'key': 'server_host',
            'value': server_host
        }
    ]
    startup_script_name = os.path.join(script_dir, 'startup.sh')
    if os.path.isfile(startup_script_name):
        startup_script = open(startup_script_name, 'r').read()
        metadata.append({
            # Startup script is automatically executed by the instance upon startup.
            'key': 'startup-script',
            'value': startup_script
        })

    shutdown_script_name = os.path.join(script_dir, 'shutdown.sh')
    if os.path.isfile(shutdown_script_name):
        shutdown_script = open(shutdown_script_name, 'r').read()
        metadata.append({
            # Startup script is automatically executed by the instance upon startup.
            'key': 'shutdown-script',
            'value': shutdown_script
        })
    config = {
        'name': name,
        'machineType': machine_type,

        # Specify the boot disk and the image to use as a source.
        'disks': [
            {
                'boot': True,
                'autoDelete': True,  # disk will be deleted together with instance
                'initializeParams': {
                    'sourceImage': source_disk_image,  # or sourceSnapshot if using snapshots
                },
                'diskSizeGb': 150,
            }
            # {
            #     'type': 'SCRATCH',
            #     'initializeParams': {
            #         'diskType': f'zones/{zone}/diskTypes/local-ssd'
            #     },
            #     'autoDelete': True,
            #     'interface': 'NVME'
            # }
        ],

        # Specify a network interface with NAT to access the public
        # internet.
        'networkInterfaces': [{
            'network': 'global/networks/default',
            'accessConfigs': [
                {'type': 'ONE_TO_ONE_NAT', 'name': 'External NAT'}
            ]
        }],

        # Allow the instance to access cloud storage, logging and all compute engine methods
        'serviceAccounts': [{
            'email': 'default',
            'scopes': [
                'https://www.googleapis.com/auth/devstorage.read_write',
                'https://www.googleapis.com/auth/logging.write'
            ]
        }],

        # Metadata is readable from the instance so we can use it to
        # pass configuration options to our instances
        'metadata': {
            'items': metadata
        }
    }
    try:
        operation = compute.instances().insert(
            project=project,
            zone=zone,
            body=config).execute()
        return operation
    except googleapiclient.errors.HttpError as err:
        if err.resp.status == 409:
            print('Instance already exists. Skipping')
        else:
            print(f'Instance couldn\'t be created: {err._get_reason()}')
    return None


def create_instances(compute, project, zone, name, count, script_dir, server_host, num_instances):
    ids = []
    for i in range(count):
        operation = send_create_request(compute, project, zone, name + '-' + str(i), script_dir, server_host,
                                        num_instances)
        if operation:
            ids.append(operation['id'])

    wait_for_all_operations(compute, project, zone, ids)


def run_command(compute, project, zone, user, name, count, command_file):
    command_content = open(command_file, 'r').read()[:-1].replace("'", "\'")
    for i in range(count):
        instance = name + "-" + str(i)
        command = ['gcloud', 'compute', 'ssh', user + '@' + instance, '--zone', zone, '--command', command_content]
        print(f'Running command\n{command}')
        out_file = open('/tmp/log-' + instance, 'w')
        time.sleep(1)  # needed because gcloud crashes miserably if run in quick succession
        subprocess.Popen(command, stdout=out_file, stderr=subprocess.STDOUT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('action', help='Action to perform: create|list|delete|stop')
    parser.add_argument('--project_id', default='metagraph', help='Google Cloud project ID.')
    parser.add_argument(
        '--zone',
        default='us-east1-b',
        help='Compute Engine zone to deploy to.')
    parser.add_argument(
        '--name', default='', help='Name (or prefix) of instances to perform the action on')
    parser.add_argument('-n', '--num_instances', default=1, type=int, choices=range(1, 200),
                        help='Number of instances to create/start')
    parser.add_argument('--script_dir', default='./',
                        help='Optional name of script to run at creation time')
    parser.add_argument('-u', '--user', default='ddanciu',
                        help='User to run commands under (for action==run)')
    parser.add_argument('--server_host', default='34.65.229.224',
                        help='The IP/hostname of the REST server that distributes jobs')

    args = parser.parse_args()

    compute = googleapiclient.discovery.build('compute', 'v1')
    if args.action == 'create' or args.action == 'c':
        create_instances(compute, args.project_id, args.zone, args.name, args.num_instances, args.script_dir,
                         args.server_host, args.num_instances)
    elif args.action == 'delete' or args.action == 'd':
        delete_instances(compute, args.project_id, args.zone, args.name, args.num_instances)
    elif args.action == 'list' or args.action == 'l':
        list_instances(compute, args.project_id, args.zone, args.name)
    elif args.action == 'stop':
        stop_instances(compute, args.project_id, args.zone, args.name, args.num_instances)
    elif args.action == 'start':
        start_instances(compute, args.project_id, args.zone, args.name, args.num_instances)
    elif args.action == 'run':
        run_command(compute, args.project_id, args.zone, args.user, args.name, args.num_instances, args.script)
    else:
        print(f'Invalid action {args.action}')
