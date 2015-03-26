#http://stackoverflow.com/questions/4395160/setting-up-github-on-a-new-computer

git workflow
=============
1. git pull
2. git add *
3. git commit -m ""
     - add commit throughout the day
     - git push throughout the day
4. for merge conflits
     - need to fix conflict and commit the results
     - make changes manually to file (HEAD: what we are on, below "=====": changes other person made)
          - delete helper things
     - git add *
     - git commit #without any options; opens multi line message
     - esc:wq 

branching
==============
1. git branch
2. git branch feature1 #completely copy of master branch
3. git checkout feature1 #switched to branch 'feature1'; make any changes and it won't affect master branch
4. make a pull request #others review 
     - one thing you want to do first is to first change to the trunk (master branch) (or whatever branch you want to merge the modified file into)
          and do a pull so the branch does not have any changes during the time that we were working on this file
          - git checkout master
          - git pull
          - if changes were made, then go back to the branch you were working on #feature1
               - git checkout feature1
          - and try to merge all changes from the trunk into your branch
               - if merge failed, then check for conflicts (edits to same file)
     - git status
     - git add *
     - git commit #Nwithout any options
     - esc:wq
     - git push

deleting file
==============
1. delete file
2. git status
3. git rm "name of file"
4. git commit -m 'message'
5. git push -u origin master

github help
==============
echo "# test" >> README.md
git init
git add README.md
git commit -m "first commit"
git remote add origin https://github.com/biohack/test.git
git push -u origin master