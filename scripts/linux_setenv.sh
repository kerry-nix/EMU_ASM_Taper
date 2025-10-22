#!/bin/sh

# Try to determine username. This helps the repo find your personal data folder.
cnt_repo_user=`id -P | sed -E 's/(.*):\*{1}(.*)/\1/'` && export cnt_repo_user
echo "export cnt_repo_user=$cnt_repo_user" >> ~/.bash_profile

# Try to determine root directory of the repo on the current file system. This helps the scripts navigate.
cnt_repo_path=`pwd | sed -E 's/(.*)\/{1}(.*)/\1/'` && export cnt_repo_path
echo "export cnt_repo_path=$cnt_repo_path" >> ~/.bash_profile

### Set up variables as temporary. This will override anything in your bash_profile.
### Use this if you want to make a temporary change, or want to avoid changing your bash_profile.
### If using this instead of bash_profile, you will need to comment out the above commands and uncomment the commands below.
### Also, you will need to run this script before using the pipeline in each terminal.

# Temporary declaration of cnt_repo_user environment variable
#cnt_repo_user=`id -P | sed -E 's/(.*):\*{1}(.*)/\1/'` && export cnt_repo_user

# Temporary declaration of cnt_repo_path environment variable
#cnt_repo_path=`pwd | sed -E 's/(.*)\/{1}(.*)/\1/'` && export cnt_repo_path

