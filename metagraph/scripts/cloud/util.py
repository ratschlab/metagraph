import os
import socket


def internal_ip():
    try:
        return socket.gethostbyname(socket.gethostname())
    except socket.gaierror:
        return '127.0.0.1'  # this usually happens on dev laptops; cloud machines work fine


def make_dir_if_needed(path):
    try:
        os.makedirs(path, exist_ok=True)
    except FileExistsError:
        pass


def dir_size_MB(dir_path):
    total_size = 0
    for dir_path, dir_names, file_names in os.walk(dir_path):
        for f in file_names:
            fp = os.path.join(dir_path, f)
            # skip if it is symbolic link
            if not os.path.islink(fp):
                total_size += os.path.getsize(fp)

    return round(total_size / 1000000, 2)


def tab(string_list):
    return '\t' + '\t'.join(string_list)


# Logger that writes data both to a file and to stdout. Can act as a drop-in replacement for a file in subprocess.Popen
class TeeLogger:
    def __init__(self, file_name, filter_str=''):
        self.log_file = open(file_name, 'w')
        self.filter_str = filter_str

    def write(self, message):
        print(message)
        if self.filter_str not in message:
            self.log_file.write(message)

    def fileno(self):
        return self.log_file.fileno()
