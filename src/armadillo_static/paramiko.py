import paramiko
ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(
            paramiko.AutoAddPolicy()) #Auto adds new hosts




#iterate over connecting with machines and passing the input command
ssh.connect(ip[iter], username='jesse', 
            password='lol')


stdin, stdout, stderr = \
        ...    ssh.exec_command("uptime")
