#!/usr/bin/env python3

import argparse
import os
import subprocess
import time

import googleapiclient.discovery


def count_pending_operations(compute, project, zone, operation_ids):
    id_filter = ''
    for id in operation_ids:
        id_filter = id_filter + '(id=' + id + ') OR '
    if len(id_filter) > 0:
        id_filter = id_filter[:-3]
    result = compute.zoneOperations().list(
        project=project,
        zone=zone,
        filter='(status!=DONE) AND ' + id_filter).execute()

    items = result.get('items')
    return items.count if items else 0


def wait_for_all_operations(compute, project, zone, ids):
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
    print('Instances in project %s and zone %s:' % (project, zone))
    for instance in instances:
        print(' - ' + instance['name'])


def send_delete_request(compute, project, zone, name):
    """ Sends a request to delete the instance with the given name """
    print('Deleting instance %s ...' % name)
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
    print('Starting instance %s ...' % name)
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
    print('Stopping instance %s ...' % name)
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


def send_create_request(compute, project, zone, name, startup_script_name):
    """ Send a request to creates a new instance with the given parameters """
    print('Creating instance %s ...' % name)
    # Get the metagraph Ubuntu 18.04 TLS image
    image_response = compute.images().getFromFamily(
        project='metagraph', family='metagraph').execute()
    source_disk_image = image_response['selfLink']

    # Configure the machine
    machine_type = "zones/%s/machineTypes/n1-standard-1" % zone
    metadata = [{
        'key': 'instance_id',
        'value': name.split('-')[-1]}]
    if startup_script_name != '':
        startup_script = open(startup_script_name, 'r').read()
        metadata.append({
            # Startup script is automatically executed by the instance upon startup.
            'key': 'startup-script',
            'value': startup_script
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
                    'sourceImage': source_disk_image,
                }
            }
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

    operation = compute.instances().insert(
        project=project,
        zone=zone,
        body=config).execute()
    return operation


def create_instances(compute, project, zone, name, count, startup_script):
    ids = []
    for i in range(count):
        operation = send_create_request(compute, project, zone, name + '-' + str(i), startup_script)
        ids.append(operation['id'])

    wait_for_all_operations(compute, project, zone, ids)


def run_command(compute, project, zone, name, count, startup_script_name):
    startup_script = open(startup_script_name, 'r').read()[:-1].replace("'", "\'")
    for i in range(count):
        instance = name + "-" + str(i)
        command = ['gcloud', 'compute',  'ssh', instance, '--zone', zone, '--command', startup_script]
        print('Running command\n%s' % command)
        out_file = open('/tmp/log-' + instance, 'w')
        subprocess.run(command, stdout=out_file, stderr=subprocess.STDOUT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('action', help='Action to perform: create|list|delete|stop')
    parser.add_argument('--project_id', default='metagraph', help='Google Cloud project ID.')
    parser.add_argument(
        '--zone',
        default='europe-west6-c',
        help='Compute Engine zone to deploy to.')
    parser.add_argument(
        '--name', default='', help='Name (or prefix) of instances to perform the action on')
    parser.add_argument('-n', '--num_instances', default=1, type=int, choices=range(1, 50),
                        help='Number of instances to create/start')
    parser.add_argument('--script', default='',
                        help='Optional name of script to run at creation time')

    args = parser.parse_args()

    compute = googleapiclient.discovery.build('compute', 'v1')
    if args.action == 'create' or args.action == 'c':
        create_instances(compute, args.project_id, args.zone, args.name, args.num_instances, args.script)
    elif args.action == 'delete' or args.action == 'd':
        delete_instances(compute, args.project_id, args.zone, args.name, args.num_instances)
    elif args.action == 'list' or args.action == 'l':
        list_instances(compute, args.project_id, args.zone, args.name)
    elif args.action == 'stop':
        stop_instances(compute, args.project_id, args.zone, args.name, args.num_instances)
    elif args.action == 'start':
        start_instances(compute, args.project_id, args.zone, args.name, args.num_instances)
    elif args.action == 'run':
        run_command(compute, args.project_id, args.zone, args.name, args.num_instances, args.script)
    else:
        print('Invalid action %s' % args.action)
