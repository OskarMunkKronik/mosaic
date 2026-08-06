#step 1: Connect to the github repository:
cd INSERT YOUR CURRENT DIRECTORY
git init
git remote add origin https://github.com/OskarMunkKronik/mosaic.git
git pull origin master

#Step 2: Install dependencies
git submodule sync --recursive
git submodule update --init --recursive
