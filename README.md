# analysis-trigger
For trigger studies

This branch is to be merged with the develop branch when running the studies with 80X samples.

After clonning the develop branch, merge this branch into develop and then reset to undo 'git add' that is done automatically with merge

cd Analysis/Trigger
git merge origin/dev-triggerstudies2017-80x --no-commit --no-ff --no-log
git reset
