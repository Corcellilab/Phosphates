import socket
import os
import sys

def find_port(timeout=2):
    """Find an available port assigned by the OS."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('', 0))  # Bind to a free port provided by the host
        
        with open('port.txt','w') as f:
            f.write(str(s.getsockname()[1]))

        return s.getsockname()[1]  # Return the assigned port number

#####

if __name__ == '__main__':
    port = find_port()

