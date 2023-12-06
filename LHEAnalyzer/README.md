## LHEAnalyzer

You need the following packages:
- [JHUGen/JHUGenMELA](https://github.com/JHUGen/JHUGenMELA) (inside the directory JHUGenMELA/ relative to the base directory)
- [MELALabs/MelaAnalytics](https://github.com/MELALabs/MelaAnalytics) (inside the directory MelaAnalytics/ relative to the base directory)
- [IvyFramework/IvyDataTools](https://github.com/IvyFramework/IvyDataTools) (inside the directory IvyFramework/IvyDataTools/ relative to the base directory)
- [IvyFramework/IvyAutoMELA](https://github.com/IvyFramework/IvyAutoMELA) (inside the directory IvyFramework/IvyAutoMELA/ relative to the base directory)

Then, please checkout this package in the base directory as
```
git clone git@github.com:usarica/LHEAnalyzer.git
```
(or with an equivalent non-SSH command).

Then, please run
```
./setup.sh -j
eval $(./setup.sh env)
```
to compile.

At the start of each session, please also run the second line above, i.e.,
```
eval $(./setup.sh env)
```
to set up the environment variables.

You can do
```
LHEAnalyzer help
```
to read the list of available options.
The executable LHEAnalyzer is located
in the executables/ directory,
which is included in the PATH environment variable after all the above.
