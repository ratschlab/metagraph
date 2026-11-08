# Cloud experiment utilities

[`remote.py`](remote.py) manages Google Compute Engine instances for the cloud
experiments. It requires Python 3, `googleapiclient`, and credentials authorized
for the target project. Review its project-specific image, network and metadata
settings before reuse. The script defaults to project `metagraph` and zone
`us-east1-b`; pass `--project_id` and `--zone` to select other targets.

The commands below describe the script's interface. Its operation polling and
`run` action have known limitations listed below.

### Usage
#### Creating instances
```
python3 ./remote.py create --name=mg-test -n=2 --script_dir=/path/to/scripts --server_host=JOB_SERVER
```
This requests `n` instances named `mg-test-0`, `mg-test-1`, and so on. Each gets
a boot disk initialized from the `mg-image-40` image in project `metagraph`.
That image project is hard-coded even when `--project_id` selects another project.

The script attaches `startup.sh` and `shutdown.sh` from `--script_dir` as instance
metadata when those files exist. Review them before creating an instance: startup
scripts run as root. Set `--server_host` to the server that distributes jobs.

`c` is accepted as an abbreviation for the `create` positional parameter, and `--num_instances` can be used instead of `-n`.

#### Starting instances
Note that creating an instance automatically starts it up, too. You only need to start an instance after explicitly stopping it.
```
python3 ./remote.py start --name=mg-test -n=2
``` 
If a startup script was attached at creation, it runs again at startup.
#### Stopping instances
Stopping an instance will shut the instance down, but its persistent data is kept (and paid for) until explicitly deleted.
```
python3 ./remote.py stop --name=mg-test -n=2
``` 

#### Deleting instances
Deleting these instances also deletes their boot disks (`autoDelete` is enabled
in `remote.py`). Copy any data you need before deleting them.
```
python3 ./remote.py delete --name=mg-test -n=2
```

`d` is accepted as an abbreviation for the `delete` positional parameter.

#### Listing instances
```
python3 ./remote.py list --name='mg-test*'
```

Lists all instances that match the given name pattern. Prefix matching is supported as in the example above.

`l` is accepted as an abbreviation for `list`.

### Known limitations

- Operation polling returns the list's `count` method instead of a numeric count
  when pending operations exist. It only checks for the absence of pending
  operations, not whether they succeeded. Check the operation result in Google
  Cloud before treating a request as successful.
- The `run` action reads `args.script`, but the argument parser does not define
  `--script`. It needs a code fix before it can execute a script. The underlying
  `run_command` function uses `gcloud compute ssh` and writes each instance's
  stdout and stderr to `/tmp/log-<instance_name>`.
