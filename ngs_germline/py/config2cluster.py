import yaml
import json
import argparse
import os

parser = argparse.ArgumentParser(description='Get the cluster configure file based the  configfile.')
parser.add_argument('--config', required=True, type=str, nargs=1,
                    help="The path of config file [required]")
parser.add_argument('--cluster', required=True, type=str, nargs=1,
                    help="The path of cluster file [required]")
args = parser.parse_args()
configfile = args.config[0]
cluster_configfile = args.cluster[0]
configfile_yaml = yaml.load(open(configfile),Loader=yaml.FullLoader)
path_log = os.path.abspath(configfile_yaml["path_log"].rstrip("/") + "/cluster_oe/").rstrip("/")
if not os.path.exists(path_log):
    os.makedirs(path_log)
# print(path_log)
cluster_config = {"__default__": {"t": "9999:00:00",
                                  "n": "1",
                                  "node": "1",
                                  "name": "{rule}.{wildcards}",
                                  "log": path_log}
                  }

for rule, thread in configfile_yaml["threads"].items():
    cluster_config[rule] = {"t": "9999:00:00",
                            "n": thread,
                            "node": 1,
                            "name": "{rule}.{wildcards}",
                            "log": path_log}
cluster_config_json = json.dumps(cluster_config, indent=4)
with open(cluster_configfile, 'w', encoding='utf-8') as f:
    f.write(cluster_config_json)

if __name__ == "__main__":
    pass
