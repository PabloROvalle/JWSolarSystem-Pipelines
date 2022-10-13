# NIRSpec pipeline installation

1. You must have the jwst package installed. An easy (and less painful) way to do that is to install previously Astroconda and then install jwst package in the same environment. JWST package, despite its hard installment in some cases, has a lot of very useful libraries for astrophysics.

2. After that, you may need to create a repo in the same directory where you have your code stored. In tha case of NIRSpec, the code will automatically create three directories stage1, stage2 and tage3 where the different pipeline levels will be stored.

3 .The variable pipestates at the begining of the code allows you to activate or deactivate some stages, if you want for example to defring you will only need the stage3 step.

4. NIRSPec pipeline is very slow, which means that during the process you will see some files created that are not in the right place. Please do not touch them and wait until the pipeline finish its process.

