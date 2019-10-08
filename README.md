# leptogenesis option

Sho's contribution to the work

> V. Brdar, A. J. Helmboldt, S. Iwamoto, K. Schmitz, *Type-I Seesaw as the Common Origin of Neutrino Mass, Baryon Asymmetry, and the Electroweak Scale*, [arXiv:1905.12634](https://arxiv.org/abs/1905.12634).

## Copyright

© Sho Iwamoto, 2019

The files under directories named `refs` are not the product of Sho Iwamoto but others, as obviously understood.
The other content in this repository is licensed by Sho Iwamoto under [the Creative Commons CC-BY-NC 4.0 International Public License](https://creativecommons.org/licenses/by-nc/4.0/).
In addition, program codes in this repository are licensed under [the MIT License](https://opensource.org/licenses/MIT).



## Usage

To install required programs, run
```sh
cd vendor
git submodule init    # initialize submodule-system.
git submodule update  # load all the submodule source codes.
./mr_build.sh         # build mr
./LoopTools_build.sh  # build LoopTools
# etc...
```

The set-up is tested in macOS with Mathematica 11.
If you encounter errors in the build process, look into each build script or consult Sho on Slack.
