import argparse
import os
import time

import googleapiclient.discovery


def list_instances(compute, project, zone):
    """ Lists the currently available instance """
    result = compute.instances().list(project=project, zone=zone).execute()
    instances = result['items'] if 'items' in result else None
    print('Instances in project %s and zone %s:' % (project, zone))
    for instance in instances:
        print(' - ' + instance['name'])


def delete_instance_no_wait(compute, project, zone, name):
    """ Deletes the instance with the given name """
    print('Deleting instance %s ...' % name)
    return compute.instances().delete(
        project=project,
        zone=zone,
        instance=name).execute()


def delete_instance(compute, project, zone, name):
    """ Deletes the instance with the given name """
    operation = delete_instance_no_wait(compute, project, zone, name)
    wait_for_operation(compute, project, zone, operation['name'])


def delete_instances(compute, project, zone, name, count):
    ids = []
    for i in range(count):
        operation = delete_instance_no_wait(compute, project, zone, name + '-' + str(i))
        ids.append(operation['id'])

    wait_for_all_operations(compute, project, zone, ids)

def start_instance(compute, project, zone, name):
    """ Starts the instance with the given name, error if instance doesn't exist """
    print('Starting instance %s ...' % name)
    operation = compute.instances().start(
        project=project,
        zone=zone,
        instance=name).execute()
    wait_for_operation(compute, project, zone, operation['name'])


def stop_instance(compute, project, zone, name):
    """ Stops the instance with the given name """
    print('Stopping instance %s ...' % name)
    operation = compute.instances().stop(
        project=project,
        zone=zone,
        instance=name).execute()
    wait_for_operation(compute, project, zone, operation['name'])


def create_instance_no_wait(compute, project, zone, name):
    """ Creates a new instance with the given parameters """
    print('Creating instance %s ...' % name)
    # Get the metagraph Ubuntu 18.04 TLS image
    image_response = compute.images().getFromFamily(
        project='metagraph', family='metagraph').execute()
    source_disk_image = image_response['selfLink']

    # Configure the machine
    machine_type = "zones/%s/machineTypes/n1-standard-1" % zone
    startup_script = open(
        os.path.join(
            os.path.dirname(__file__), 'startup_script.sh'), 'r').read()

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

        # Allow the instance to access cloud storage and logging.
        'serviceAccounts': [{
            'email': 'default',
            'scopes': [
                'https://www.googleapis.com/auth/devstorage.read_write',
                'https://www.googleapis.com/auth/logging.write'
            ]
        }],

        # Metadata is readable from the instance and allows you to
        # pass configuration from deployment scripts to instances.
        'metadata': {
            'items': [{
                # Startup script is automatically executed by the
                # instance upon startup.
                'key': 'startup-script',
                'value': startup_script
            }]
        }
    }

    operation = compute.instances().insert(
        project=project,
        zone=zone,
        body=config).execute()
    return operation


def create_instance(compute, project, zone, name):
    operation = create_instance_no_wait(compute, project, zone, name)
    wait_for_operation(compute, project, zone, operation['name'])


def create_instances(compute, project, zone, name, count):
    ids = []
    for i in range(count):
        operation = create_instance_no_wait(compute, project, zone, name + '-' + str(i))
        ids.append(operation['id'])

    wait_for_all_operations(compute, project, zone, ids)


def wait_for_all_operations(compute, project, zone, ids):
    print('Waiting for *all* operations to finish...', end='', flush=True)
    sleep_time_sec = 1
    while True:
        if count_pending_operations(compute, project, zone, ids) == 0:
            print("done.")
            return
        print('.', end='', flush=True)
        time.sleep(sleep_time_sec)
        if sleep_time_sec<16:
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
        '--name', default='metagraph', help='Name (or prefix) of instances to perform the action on')
    parser.add_argument('-n', '--num_instances', default=1, type=int, choices=range(1, 50),
                        help='Number of instances to create/start')

    args = parser.parse_args()

    compute = googleapiclient.discovery.build('compute', 'v1')
    if args.action == 'create' or args.action == 'c':
        create_instances(compute, args.project_id, args.zone, args.name, args.num_instances)
    elif args.action == 'delete' or args.action == 'd':
        delete_instances(compute, args.project_id, args.zone, args.name, args.num_instances)
    elif args.action == 'list' or args.action == 'l':
        list_instances(compute, args.project_id, args.zone)
    elif args.action == 'stop':
        stop_instance(compute, args.project_id, args.zone, args.name)
    elif args.action == 'start':
        start_instance(compute, args.project_id, args.zone, args.name)
    else:
        print('Invalid action %s' % args.action)
