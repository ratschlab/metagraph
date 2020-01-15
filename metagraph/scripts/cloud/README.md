# Getting ready
### Install the Google Cloud SDK
  - Install python3 (I used 3.7)
  - Install the [Google Cloud SDK](https://cloud.google.com/sdk/docs/quickstart-macos) as described here, skip authenticating.
  - Create a [service account](https://cloud.google.com/docs/authentication/getting-started) for authentication

At installation of the Google Cloud SDK, When asked to set a default compute region, enter `europe-west1-c`. This is where
our instances are located now (in Zurich!).

### Usage
#### Creating instances
```
python ./main.py create --name=mg-test -n=2 --script=<startup_script>
```
This creates `n` instances (running Ubuntu 18.04) with the persistent disk `metagraph-sample`, which contains all the code, 
binaries and deps needed to run `metagraph`. 

At startup, each instance will run the shell script optionally passed in via the `--script` parameter (as root!).

`c` is accepted as an abbreviation for the `create` positional parameter, and `--num_instances` can be used instead of `-n`.

#### Starting instances
Note that creating an instance automatically starts it up, too. You only need to start an instance after explicitly stopping it.
```
python ./main.py start --name=mg-test -n=2
``` 
If the a startup script was specified at creation using the `--script` argument, this script will be executed again at startup.
#### Stopping instances
Stopping an instance will shut the instance down, but its persistent data is kept (and paid for) until explicitly deleted.
```
python ./main.py stop --name=mg-test -n=2
``` 

#### Deleting instances
Use this when you are done with the instances and no data produced by the instances on the local hard drives is needed. 
Deleted instances are lost forever.
```
python ./main.py delete --name=mg-test -n=2
```

`d` is accepted as an abbreviation for the `delete` positional parameter.

#### Listing instances
```
python ./main.py list --name=mg-test*
```

Lists all instances that match the given name pattern. Prefix matching is supported as in the example above.

`l` is accepted an abbreviation for `list`.

#### Executing a script
```
python ./main.py run --script=[path_to_script] --name=mg-test -n=2
```

Runs the given script on n instances with the given name prefix. The output (stderr and stdout) of the executed script
for each machine is written locally to `/tmp/log-<instance_name>`. 