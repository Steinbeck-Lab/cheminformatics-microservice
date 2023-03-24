# Changelog

## [0.7.0](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/compare/v0.6.0...v0.7.0) (2023-03-24)


### Features

* add CDK aromatic ring calculator [#84](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/issues/84) ([0f762ae](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/0f762ae355b54d93228f4e052b7cdd2063a8dfba))


### Bug Fixes

* [#81](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/issues/81) POST requests now parse body using fastapi interface and swagger UI ([47c9908](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/47c9908c23575ab105f603fdb46b0dfa9bc70886))
* docker build command in prod ([c3433f9](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/c3433f94e2e2ee570de59d8713a5d8646ae9dccf))
* move test job before creating a release. ([d4679cc](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/d4679cc5b4fd308db53735d725b5baf1e762da3e)), closes [#88](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/issues/88)
* updated default release version from "pre-release" to "latest" ([e086b19](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/e086b19f478ec7a64f57bcf15b57bbddeb46b736))

## [0.6.0](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/compare/v0.5.3...v0.6.0) (2023-03-24)


### Features

* add All descriptor module for CDK and RDKit [#79](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/issues/79) ([0e28187](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/0e28187ae6a3e64eb1ed8de12aa80d496839058b))
* add CDK descriptors for COCONUT [#78](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/issues/78) ([d850fde](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/d850fde603960d22cd1637cc5c8a0e2c40cefd9e))
* add function to COCONUT descriptors based on toolkit [#80](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/issues/80) ([a6db174](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/a6db1747bc712c1b6ff317b4dab12b2f3fe3461f))
* Add table preview & JSON,HTML toggle [#76](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/issues/76),[#77](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/issues/77) ([7e2e063](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/7e2e0632f22f1bccee02dcbed387a6c67db8b3cc))


### Bug Fixes

* add nplikeliness for descriptors [#74](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/issues/74) ([41eaa94](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/41eaa94436db835bd94918e55cdeb67c50e09995))
* descriptor dictionary order ([fcdb1cb](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/fcdb1cb622fd5aab162b8fc60ef3811b846b75d2))
* remove the if cond to check branch name ([85fcda0](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/85fcda039ca34c7b925c31c011af446049b9f924))

## [0.5.3](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/compare/v0.5.2...v0.5.3) (2023-03-23)


### Bug Fixes

* remove other release cond ([9e21ac0](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/9e21ac027f410d3aa32123aa1c64759be24ad730))

## [0.5.2](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/compare/v0.5.1...v0.5.2) (2023-03-22)


### Bug Fixes

* add CDK badge ([53f919a](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/53f919a58278c8dcc4085ab3f389f9a10d7d3060))
* update token in release-please to use PAT ([12720ce](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/12720cea837c68c6f88dea3499d9e989c9a03906))

## [0.5.1](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/compare/v0.5.0...v0.5.1) (2023-03-22)


### Bug Fixes

* add prerelease condition to prod build ([fb34d7f](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/fb34d7fc1151e3e741356238bd775657c487ebc9))

## [0.5.0](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/compare/v0.4.0...v0.5.0) (2023-03-22)


### Features

* add SRU, Murko framework and implement COCONUT related descriptors ([abf9991](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/abf9991815adbc98c48a0c1f19fdf1e0d162d326))


### Bug Fixes

* docker build command ([3da9cfb](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/3da9cfb20f39ba189049dc5009e70e550b77dbbe))
* docker build syntax ([5b3f50f](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/5b3f50f7fa76ca6bceeee9107121b453f0f8e419))
* dynamic version number fed into the swagger api from docker release ([961aa11](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/961aa119c54cc3e58c5c9b6d932b94fcfbccd9df))
* fetch pre-release as well ([71d35fd](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/71d35fd4635a79e0481cda144a873b24ff33c03b))
* fetch release tag while build ([0727dfc](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/0727dfc02f9dc30d3f594a078e12a7a48bbca228))
* indentation issues ([7cf82d4](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/7cf82d477311e4179c1c380ac6c36d9d84ad1cf1))
* linting issue resolution ([8537b83](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/8537b83d30084dff5c577c7ec86dbf92f3d86885))
* module import ([fa7f5a8](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/fa7f5a8c38e456d8b4e5ee7340488feb5112bb72))
* remove unnecessary instalaltion in dev-build ([245a188](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/245a1887bdf211b1f390803f4d7ec13fc89cf544))
* SRU import ([988f76c](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/988f76c83c58d3ccb035c251544e0945bb14b029))
* update trigger event for prod build workflow. ([e8cde48](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/e8cde48beb5122765a76dfbd606b83fb80bc6a71))

## [0.4.0](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/compare/v0.3.0...v0.4.0) (2023-03-21)


### Features

* Add 3D depiction module ([5767a49](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/5767a495e6c0522d0c1e68a7ff1cbbf7358d4676))
* add 3D visualization command to ReadMe ([7791a81](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/7791a81afcb29ddf23020240650642e0d5417d71))


### Bug Fixes

* add Zenodo entry and citation ([aa4f733](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/aa4f7334ea2380489dd1dae40b3866ed36531a5c))
* enabled 3d viewer ([043d4ad](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/043d4ad658409703049ac05f05fc6f56fb7c2947))
* linted code ([1e63fa4](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/1e63fa4453b0a5f6175e424c298bc5ac710e12fe))
* Linter issues on depict3D ([0626491](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/0626491680275ed0a82edb5503f01dccf1445554))
* linting error fix ([185ec4c](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/185ec4c6089d0a0bc9a632a3daf06ba686ff4aef))
* logo padded, resolves [#52](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/issues/52) ([df21d02](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/df21d0230dd3d7389374982b8573f2572c116344))
* remove py3Dmol from requirements ([6ca6931](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/6ca69318116efdeddeaa87d5fa45a380822e93c1))
* remove unused modules and imports ([1b7e1fe](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/1b7e1fe2352bcec88e73e1216f96efdc5dbc4f1e))
* Removed matrix testing on prod and dev. Removed tests on dev build ([7268750](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/7268750b445b183e4ba05f67b6a2489a5febb43f))
* update release version ([7945fbb](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/7945fbb845a1fe7b8363ef428870defc59f580bb))
* updated readme urls to point to production, redirected root to docs ([7e2a5ae](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/7e2a5ae8da3ba965441ea6d383aa4c2c752a4de5))

## [0.3.0](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/compare/v0.2.2...v0.3.0) (2023-03-18)


### Features

* customised fastapi swagger interface ([0b1f620](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/0b1f620296cff5b1b76fff291492650fd518f912))


### Bug Fixes

* scaling and atom colouring ([c6440bd](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/c6440bda7deef608d16e9c27b719900ebc69de6b))

## [0.2.2](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/compare/v0.2.1...v0.2.2) (2023-03-17)


### Bug Fixes

* links to examples ([73144f0](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/73144f07b4ee3cb5d7d451bb7821744b5de72856))

## [0.2.1](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/compare/v0.2.0...v0.2.1) (2023-03-17)


### Bug Fixes

* added Decimer routes to the main ([b7df835](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/b7df835707d3424923f91c5748be1310e67222ce))

## [0.2.0](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/compare/v0.1.0...v0.2.0) (2023-03-17)


### Features

* add converters and decimer modules ([91778b5](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/91778b5e4639e1f15e5da91df183cd9374405060))
* add documentation ([fab28a9](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/fab28a98d75aa7d9c37f9763738d35f82255887c))
* add logo ([f24a7ce](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/f24a7ce042789aa7b994a81968e3a83e164a8f80))
* add MMFF94 optimization ([4ca0d01](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/4ca0d01e225da6868c9714d1448c0f8857f82372))
* add pytest ([259f404](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/259f4046e9e2f5771ae70469ebd5e7ada421d4ab))
* add pytest to dev-build ([ae067d8](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/ae067d8acb9a64fae4da6574794e5243ec32f2d8))
* add RDKit random 3D conformer generation ([c78553b](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/c78553b7aa17d4d7a0065ccef4698452b5276c22))
* add ReadMe ([d426803](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/d42680393258a4e94296f88e575473f2362dc5d1))
* added code comment blocks, added additional endpoints for inchi, inchikey and cannonicalisation, merges smiles conversion from iupac / selfies ([81a6d3e](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/81a6d3e26090e544adbe3f866fe55cb619d18540))
* Create CITATION.cff ([de4289d](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/de4289d35b4915685a3afeab51ce09155899f084))
* Initial tests and fix selfies decoder ([839e908](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/839e9084564039bad4893386d61f578a096785cc))
* linted code, merged mol end point and added extra query parameter to choose generator ([32cc977](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/32cc9779641a28fe2ecde039e780b137a59fe9c2))


### Bug Fixes

* add citation ([d0ec378](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/d0ec3781c4fee65549bbef79e2fc85534709d229))
* file name mismatch issue fix (502 error) ([9fe53e5](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/9fe53e5ce39d0f2d550ac409a4ca9395c0c402d5))
* GIF image conversion issue [#38](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/issues/38) ([51817b8](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/51817b8ce0c41e56ebca526d5641b08a506c8f65))
* hyperlinks of university ([007631b](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/007631b6ff0bd80459f5d9095eac27964af094c8))
* install pytest ([a2a257b](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/a2a257b72bb6a8aa0ff963cd5f85b40287613232))
* linter errors ([b720837](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/b7208375b24190ccd1150a84523bafaa1a80fe2a))
* remove database.py ([dc59266](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/dc59266c057ee5b0265667cb03486b89670fca6f))
* remove requirement installtion ([cf34a2c](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/cf34a2cdccda613c30ab09345141858bee19a701))
* rename converters to convert ([563ceb7](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/563ceb7bf15b5d9f84e3426fc029ce1dea75027d))
* replace branch reference ([04fbc92](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/04fbc92596140e7b29c0407a39bf8b07fcbc0a10))
* resolve pytest issue ([2135952](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/21359526880de508866b27cdc2ead6bd815a3e37))
* resolved import convert error ([421847b](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/421847b711a95a7e942beee29a69b81d862c685a))
* Test indentation and linting issues ([424e09d](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/424e09d093c77ce57b19c8c79e5a0f2f72e501d4))
* test installation ([50f10f0](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/50f10f0247cef513bafd43ee1561b75a62d90f27))

## 0.1.0 (2023-03-14)


### Features

* add build for prod and seperate docker file ([d85d726](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/d85d726539b59b84ddcadbe3257faefc56a9859f))
* add docker-compose file ([103df63](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/103df63027ce4a05f62f0103ca624eaed17b251f))
* add github workflow to ssh into remote ([2c0d3cc](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/2c0d3cc2c0b71b67eb29c60d10601a1bc80bf7d1))
* add jobs to build,test & publish. ([e9ac425](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/e9ac425d042f1fbdba26d81a203b135779d004e5))
* add release-please ([7eda0ca](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/7eda0ca4d5abe117ccd3b5b31086b122b4b23862))
* add roatation to CDKDepict ([01b90f0](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/01b90f06efa35459621c593f88ad18d130a0cb4f))
* add rotation to CDKDepict ([9e3c29e](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/9e3c29ed6de98606d74ac74b3274721c3bf946d9))
* add test workflow ([d484b9d](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/d484b9d347b967883cc9c121f4ced6f6e7b63501))
* added classyfire module ([76295c7](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/76295c7036e701be471a602cddc41bad5d717149))
* added decimer image processing end points ([6ccf663](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/6ccf663f618ae3ae79bfb50e2794f3fe0e72f5a2))
* added endpoints to expose cdk and rdkit depiction routines ([f8a05ba](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/f8a05ba5afeb3afdc402719d26945aad1da0355e))
* added iupac to smiles conversion end point ([eb17cae](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/eb17cae94fbdb60078fa375d0b92aea37de77456))
* added latest RDKit version, psycopg2, rdkit cartridge and various other updates ([746baaf](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/746baaf63d93519c4eadfda680924184e01c36f9))
* Added support for CDK SDG 2D coordinates ([b3acfe1](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/b3acfe16edb5d2d77e0c45a67bc083e150d8f769))
* build and publish docker image ([7a4aa0e](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/7a4aa0eda2bb85de7bcdc335172782a29d75ae28))
* CDK and RDKit depict functions ([92c269a](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/92c269a8ddba49f38d1378c20aba5632f73df673))
* deployed stout package ([416fb7b](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/416fb7ba8d428a09ea24dd8c42d8a5c3fd4e8665))
* descriptor calculations using RDKit ([7b836ed](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/7b836ed4d5bedbc478b6fd91f880791ff488bfc0))
* implement RdKit based np likeliness score ([2c5a639](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/2c5a639703d15a2c2c95390f3909265036ad4858))
* selfies encoding and decoding ([98b6547](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/98b6547c2c1610db10b21f534a621e0a7cf63b81))
* update requirements, stop installing dependencies twice and added decimer segmentation ([096dc67](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/096dc671f7326a4be85d541970ce75175d92fc2c))


### Bug Fixes

* "Unsupported upgrade request issue: uvicorn" bug fix ([fe5fa28](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/fe5fa28ea0f69625ef9cea4278bb0235ef2bf71a))
* add rotation to CDKDepict ([be6cb4a](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/be6cb4a1bb71606081277b26ca5386e56d813d9d))
* added websockets to requirements doc to resolve "Unsupported upgrade request." ([0e3a0f6](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/0e3a0f63793955463c17d6137d716270508388d6))
* added yes flag for confirmation prompt ([2214bdf](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/2214bdf76a026307300dc425932e59548ed1402d))
* all flake8 errors ([81f2e05](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/81f2e05a4fdaca067bafe36098cdecfb51d9deda))
* bumped stout to latest version ([dde0ad8](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/dde0ad80acbc09e8c14d457f1cabf9eb5ee973c6))
* CDK2D string display ([86938fa](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/86938fa5c7f350dc86f8391346a67abbbffd0b7e))
* Docker dependencies ([c53917c](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/c53917ca054a33f7e59d2836d434bf543dc315af))
* enabled incoming image via post request ([4e168cb](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/4e168cbe99d6405e32ff451debfff77bd26c22fc))
* enabled null img uri exception handling on server end ([c49e76f](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/c49e76fd0c7065b21024d376adc0d48387f04bf2))
* exposed geometric rotation parameter via query param ([0f481f9](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/0f481f9f6855889c332610f3804e4bb2794a4aa4))
* flake8 issues ([fd596da](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/fd596daa483aad3337ddcaa45e9c50d4089a3083))
* flake8 linting issues ([63d03b7](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/63d03b72ba546ec0884f8a213e84c2b39930e87c))
* Import np score function ([4b8b351](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/4b8b3512c2961f83c84f29ebb72e0c7cd8aa38c1))
* re-add the build workflow ([024ff15](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/024ff1542eac3cd32e87e7a23dccc74725cadd95))
* refer to development branch ([8ad9031](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/8ad9031847a44c919b1f9df0f68cb162f021938f))
* remove build.yml ([8eb5f4b](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/8eb5f4b92adb6a3d294cd185bec035039bb8226e))
* remove build.yml file ([3cd0914](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/3cd09146eefaf6ab089130d152efaf1b12857e66))
* remove Cluster configuration ([46e7f81](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/46e7f81fe6caaa65da7848ee17f67573cf67e491))
* remove RDKit build install via PyPi ([c484b7c](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/c484b7c15b2a347e47a51c7e9d59eb196370712a))
* removed shiftdb import scripts ([3eb90d3](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/3eb90d3cd1bfb853ab6c7760adfce9a8c45610fc))
* resolve conflict ([560287a](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/560287ac79a2d59d812dfa61cd97532a90d9b3cc))
* saving file to tmp folder and also using original file name ([cb73dcb](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/cb73dcb20c33bfb1fb429fdeecb36e136b0ee1e7))
* Scaling CDK SVG Depict ([5b6f2fc](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/5b6f2fcd0fd805aace8f8caf73d237b6d8968899))
* try to create file in rmeote server ([f809971](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/f8099711e1d82f1c4b7df5a2757280f859953f2b))
* update docker repo name for dev ([594f3c9](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/594f3c9a6167468edc78a6060318cfef88a39c06))
* update login to remote server ([1e4e6ad](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/1e4e6ad25b82d75b24db69f1f91d9db7a38ae61c))
* update python version ([1031404](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/103140492189cac53b2969b020c7fb90381db46a))
* update repo name ([37657e0](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/37657e0ebc0e027d67d5300c9523fec72874ce96))
* update the if condition ([707848f](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/commit/707848fe4e2e4d9fcb4e92f38a3227b442a62714))
