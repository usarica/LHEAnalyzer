LHEAnalyzer

You need the following packages:
- JHUGen/JHUGenMELA
- MELALabs/MelaAnalytics
- IvyFramework/IvyDataTools
- IvyFramework/IvyAutoMELA

Then, do
```
./setup.sh -j
eval $(./setup.sh env)
```
to compile.

At each session, please also run the second line above, i.e.,
```
eval $(./setup.sh env)
```
to set up environment variables.

You can do
```
LHEAnalyzer help
```
to read the list of available options.
The executable LHEAnalyzer is located
in the executables/ directory,
which is included in the PATH environment variable after all the above.
