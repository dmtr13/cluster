#!bin/bash
cd .. && bash push.sh
read -p "Commit message: " desc
cd wiki
git add .
git commit -m "$desc"
git push
cd ..
git add .
git commit -m "$desc"
git push
cd bin
