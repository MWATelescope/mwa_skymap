# mwa_skymap
All-sky maps and movies showing the MWA telescope primary beam

## Publishing and Releasing

There is a github workflow CI called release.yml which will:
* Run tests
* Build the package
* Publish to PyPi (Using the "test" or "real" PyPi- currently hardcoded to "test")
* Create a GitHub release

Before first use you need to:
* In the GitHub repo, create a new environment called "pypi". 
  * Navigate to your repo on GitHub → Settings → Environments (in the left sidebar under "Code and automation").
  * Click New environment, type pypi as the name, and click Configure environment.
  * Optionall add protection rules- 
    * Required reviewers — adds a manual approval step before the publish job runs. Good for preventing an accidental tag push from immediately shipping to PyPI.
    * Wait timer — adds a delay before the job proceeds.
    * Deployment branches and tags — restrict deployments to only tags matching v* so the environment can't be used from a random branch.