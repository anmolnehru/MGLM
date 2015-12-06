#! /bin/sh
#Author - Anmol Mohanty

#creates key file called brain_aws.pem for ssh sessions
aws ec2 create-key-pair --key-name brain_aws --query 'KeyMaterial' --output text > brain_aws.pem

#changes permissions || Imp step, don't remove
chmod 400 brain_aws.pem

#spawns 19 instances(limit of 20, inclusive) of type c4.8xlarge with security brain_aws file & pipes output to servers.json
aws ec2 run-instances --image-id ami-60b6c60a --security-group-ids sg-35355d53  --count 19 --instance-type t2.large --key-name brain_aws


#IMPLEMENT A PAUSE TO WAIT UNTIL INSTANCES ARE BOOTED UP
#sys.pause(10);
sleep 2m #wait 2 mins


#describe the instances and pull the public ips(non elastic) || maybe better to do in prev step only??
aws ec2 describe-instances > servers.txt
grep "PublicIpAddress" servers.txt > public_ip.txt
#should store list of ips in public_ip.txt ##done


#interesting nice script to properly format the file
vim -S vim.trim public_ip.txt

#remove the master ip address from this list
sed '/52.70.48.243/d' filename.txt


#parse servers.json? (python)? to obtain list of private ips of children
#ip_list=parse(servers.json) #python call


#STATE RUNNING
#servers running state

#iter=0 #iterates over servers
#
#while [$a -lt 19]
#do
#	#secure copy of data partition from local to remote
#	scp -i brain_aws.pem data[$a] ip_list[$a]:~/
#	#explore possibilities with paramiko, workload distribution
#	#python potentially easier
#done

#start the 20th job in the backgroud
./shared/legr_dti_parallel 20/. shared/. &

#pulls up the list of ips and connects to them seamlessly
#and copies over various portions of data
./execute_command_from_ipfile.sh public_ip.txt
#better to be in blocking mode for safety, anyway not too long

#executes individual queries (p_value computations) on slaves parallely
./parallel_execute_slaves.sh public_ip.txt

#move the local file to the final resting place
mv p_value.txt ~/PVALUE/20

#Upon this job completion messages will come in and all pvalues should get accumulated
#begin verficiation of p_values

#verify by number of values and value < 1, can be done in the c++ file too; think about this
verify(p_values)

#will be interesting to see the various cases
if(failure in an index) -- respawn that directory on a new instance



#job mostly copmpleted, will wait for the p_value texts to come flooding in
#stitch them together?




