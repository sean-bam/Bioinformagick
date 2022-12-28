layout: page
title: "Version control with Git"
permalink: /Git

## <a name='TOC'>Table of Contents</a>

1. [Setup a Git repository from the cmd line](#Setup)
2. [Revert a specific file to a commit](#Revert)


**<a name="Setup">Setup a Git repository from the cmd line</a>**

Useful [instructions](https://docs.github.com/en/github/importing-your-projects-to-github/adding-an-existing-project-to-github-using-the-command-line) to setup a repo from the command line
```bash
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/sean-bam/**REPO**.git
git push -u origin main
```

**<a name="Revert">Revert a specific file to a commit</a>**

Revert a specific file to a previous commit `c5f567` 

[Source](https://stackoverflow.com/questions/215718/how-can-i-reset-or-revert-a-file-to-a-specific-revision)
```bash
git checkout c5f567 -- file1/to/restore
```
