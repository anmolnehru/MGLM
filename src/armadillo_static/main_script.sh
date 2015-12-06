#! /bin/sh
#Author - Anmol Mohanty

#creates key file called brain_aws.pem for ssh sessions
aws ec2 create-key-pair --key-name brain_aws --query 'KeyMaterial' --output text > brain_aws.pem

#changes permissions || Imp step, don't remove
chmod 400 brain_aws.pem

#spawns 19 instances(limit of 20, inclusive) of type c4.8xlarge with security brain_aws file & pipes output to servers.json
aws ec2 run-instances --image-id ami-60b6c60a --security-group-ids sg-35355d53  --count 20 --instance-type c4.8xlarge --key-name brain_aws
#describe the instances and pull the public ips(non elastic) || maybe better to do in prev step only??
aws ec2 describe-instances > server.txt
grep "PublicIpAddress" servers.txt > public_ip.txt
#should store list of ips in public_ip.txt ##done
vim -S vim.trim public_ip.txt


#parse servers.json? (python)? to obtain list of private ips of children
ip_list=parse(servers.json) #python call

#takes few minutes for servers to spin up
#implement a sleep or pause mechanism



#IMPLEMENT A PAUSE TO WAIT UNTIL INSTANCES ARE BOOTED UP

#STATE RUNNING
#servers running state

iter=0 #iterates over servers

while [$a -lt 19]
do
	#secure copy of data partition from local to remote
	scp -i brain_aws.pem data[$a] ip_list[$a]:~/
	#explore possibilities with paramiko, workload distribution
	#python potentially easier
done

#job mostly copmpleted, will wait for the p_value texts to come flooding in


#stitch them together?




