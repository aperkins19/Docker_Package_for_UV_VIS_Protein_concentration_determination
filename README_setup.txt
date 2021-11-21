This generic app will allow you to easily export and run your data analysis scripts on remote machines.

SETUP

1. Define in the dockerfile which python version you want to use.
2. Define which python packages your script needs in requirements.txt
Be careful with package dependancies - it's best to define the version of each package.

3. Name your script main_script.py

LOAD THE APP ONTO THE REMOTE MACHINE.

1. Ensure you are connected to the university VPN.
2. Ping the remote machine to ensure you are in the system properly.

3. If it is your first time working on the machine, SSH into the machine using your university username to establish the correct keys
e.g.

ssh s1530400@SCE-BIO-C06156

4. Logout and back in using sftp

sftp s1530400@SCE-BIO-C06156

5. Navigate both to the local and remote directories of where the app is and where you want it to go.
e.g.:

Remote -> cd desktop
Local  -> lcd desktop

6. Make a directory on the remote computer to put the app into.
e.g.

mkdir remote_number_crunching

7. Navigate inside
e.g.

cd remote_number_crunching

8. Upload the app by transfering all the files ("*") in the local directory.
e.g.

put *

9. Check it's there.
e.g.

ls


BUILD AND RUN DOCKER CONTAINER.

1. Login using SSH

2. Check Docker is installed.
e.g.

docker --version

3. Navigate to the folder remote_number_crunching

4. Build the docker image. This will take a long time if you have lots of dependancies and packages.
e.g.

docker build -t my-python-app .

5. Run the docker image.
e.g.

docker run -p 8000:5000 -it my-python-app



https://www.youtube.com/watch?v=YFUhdxI4kcA