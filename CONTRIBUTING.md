# Workflow

Below is an outline of our general workflow. It is primarily intended for our software team, 
but outside contributors may wish to draw from it as well.

1. Using gitflow locally, start a branch and make any desired changes.
2. Update the RELEASENOTES.md document and push the branch to GitHub.
3. Create a PR that has a base branch of develop (for features and bugfixes) or master (for releases and hotfixes).
4. If the checks pass and we're happy with the release, merge the PR. **Do not delete the branch on GitHub yet!**
5. On the local copy, run `git flow finish` for the branch you started. (At this point, it is safe to delete the branch on GitHub.)
    * If you merged back to develop, this is the last step. If you merged back to master:
6. Push the local copy of develop.
7. Push tags.
8. Create a release on GitHub from the version tag you just added.

# Standards

All features of `galactic-spin-W1` should be tested with Travis-CI wherever possible.  This is done
through the top-level `.travis.yml` file.  You can find Travis-CI's docs
[here](https://docs.travis-ci.com/).

Addition of features should be accompanied by some form of documentation outlining the use of the
new features. Similarly, if a feature is altered, documentation should be adjusted to match.

Variable and function names should relate to what the variable or function is doing.  
**Unacceptable:** `func` - Name does not describe what the function does at all.  
**Acceptable:** `calculate` - The name indicates a category (that the function calculates something), but it is still vague.  
**Best:** `integrate` - The name clearly indicates what the function does (take an integral).

*(Spot something that's not up to our standards? Submit an issue or pull request!)*

## Python Standards

Keep code neatly organized:
* Leave space before and after multi-line function definitions
* For jupyter notebooks, try to break up cells where it is sensible - for example, separating plotting cells, cells that define functions and variables, and cells that process information.
        * In general, avoid having notebooks that consist only of a few very long cells or entirely of very many very short cells.
* Avoid long lines of code where it's possible to break the code up into multiple lines instead.

Example:

```python
def f(a, b):
    return a*b

i = 0
if (1 == 1):
    while (n < 3):
        i += 1
        f(1,1)
```

## Documentation Standards

All functions accessible to the user should have at least one example provided, possibly more if
the usage is complex or varies significantly.

All variables and options available to the user should be clearly defined.

The purpose of each notebook should be clearly defined somewhere. 
If notebooks are intended to be run in a specific order, that information should be documented as well,
and the notebooks should have a numbering system to indicate this.

Documentation files should, at the end of the file, note the date corresponding to the last time
they were updated as well as the relevant version number.


# Pull Requests

Pull Requests (PRs) are created in order to submit to the owner(s) of the repository some code for
consideration. Typically that code will address some issue or improve the code in some way, we
should be clear about how we expect PRs to improve the code in our contributing documentation.
When creating the pull request you have to supply a comparison branch.  When submitting PRs,
please default to submitting to the `develop` branch and never submit directly to the `master`
branch.  This allows us to correctly manage versioning.  If you submit a pull request directly to
the `master` branch, you will be asked to change the target to `develop` (or another applicable
branch).

In order for us to be consistent with GitFlow we should adhere to the following:

**PRs submitted by the development team:** we should locally create a feature branch using the
gitflow protocol on our local machine, and then push that branch. That branch should then be
selected as the "comparison" branch for the PR. Further for the merger to be compatible with
gitflow we should define the base branch as "develop." So the steps are:

1. create a feature branch locally make some changes and push it to the remote (GitHub)
2. open a pull request with base branch `develop` and comparison branch the feature you just created `feature/XXX`
3. When you're done committing (and pushing!) to the feature branch push the button on GitHub to merge the PR back--it will merge it to develop 
4. Delete the branch on github and your local machine and add notes to the upcoming release
5. the feature will be released when the code team does the next release. 

**PR requests submitted from outside our development team:** are very similar to those from the
development team, but the team won't have access to or control over the feature branch created. It
would be created by a fork of the repository. So it looks like this:

1. create a fork of the repository with a branch dedicated to the issue (could be the local `master` we can't enforce any naming conventions there). 
2. open a PR with base branch `develop` and the comparison branch the branch on the fork you just created. 
3. When you're done committing alert the development team in the PR by using the @villaa or other tags. 
4. This will be merged back by the development team if the criteria for code improvement are met. 
5. the feature will be released when the code team does the next release. 

All PRs will be automatically by Travis-CI.  Please note whether you updated the CI or
whether no change was needed.  If for some reason a new, untested feature is implemented, but you
are unable to implement the necessary CI, explain why and how it can be manually tested.

## Release Documentation

When a PR is accepted it will be staged for a release. We will make note of all the currently
staged changes in the RELEASENOTES.md file. It is helpful, but not necessary to put a short
description under the `Next Release` section briefly describing the work done in a PR.  
A template is provided in `pull_request_template.md`.

**Other information:**
Anything else you want to say.

# GitHub Issues

Issues fall into three categories:
* Bug report
* Feature request
* Documentation issue

When submitting issues, please be specific. We can only resolve a bug if we know what about the
program isn't working, implement a feature if we know what aspect is being improved, and clarify
documentation if we know what part is unclear.

Below are outlines for determining what your issue qualifies as. When
submitting an issue, please specify which of these three categories you think it belongs in. We
understand that the three categories can overlap, so don't worry too much if you aren't sure if
the category you chose is appropriate. When creating an issue, you will be given the option to 
select a template; these are just to
help people know what to write, and their use is not strictly required (although it may help us
address the issue faster).

## Bug report

When submitting a bug report, please make sure to include any information necessary for us to
reproduce the bug. If we can't reproduce it, it will be much harder to diagnose and solve the
issue.

An issue is a bug report if:
* The code does not run or only partially runs.
* The code does not build.
* A command associated with the code fails to run despite matching the documentation.
* The code takes an inordinately long amount of time to run, taking hardware into account.
* The code gives no output either to a binary file, a log file, or the terminal, when it should be giving some kind of output there.
* A command that should give consistent output gives different output each time.
* The result of a command directly contradicts what the documentation says should occur.

An issue is not a bug report if:
* The code does not interface with an environment that the documentation does not specify it will interface with. (Feature request)
* The code is missing the ability to do something you think it should be able to do, but the documentation does not specify it is capable of. (Feature request)
* The documentation is unclear but the code does not give results that directly contradict it. (Documentation issue)

## Feature request

An issue is a feature request if:
* You are requesting for the code to interface in a new way that it currently does not, such as a new command or argument.
* You are proposing a particular way to increase the speed of the code.
* You are pointing out where the code could be more user-friendly.
* You are otherwise requesting for the code to do something it is not yet written to do.

An issue is not a feature request if:
* It does not affect the code, only the documentation. (Documentation issue)
* It is to fix unexpected behavior. (Bug report)
* You are providing the feature you are requesting. (Pull request)


## Documentation

An issue is a documentation issue if:
* A command mentioned in the documentation cannot be found.
* You would like an example and there is no similar example in the documentation.
* There is a part of the documentation you are asking for us to clarify.
* There are spelling and grammar errors, missing images, or broken links which you do not know how best to fix with a pull request.

An issue is not a documentation issue if:
* You provide fixed wording, spelling, and/or grammar for all issues you point out. (Pull request)
* The code attempts to run but fails. (Bug report)
* You are looking for a way to do something that you do not know exists and is not mentioned in the documentation. (Feature request)

*Last updated 21 January, 2022, v0.1.0*
