#step 1: Connect to the github repository:
cd INSERT YOUR CURRENT DIRECTORY
git init
git remote add origin https://github.com/OskarMunkKronik/mosaic.git
git pull origin master

#Step 2: Install dependencies
git submodule add https://github.com/OskarMunkKronik/regionofinterest.git
git submodule add https://gitlab.com/tensors/tensor_toolbox.git
git pull --recurse-submodules



