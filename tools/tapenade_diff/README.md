# tapenade diff

Tapenade is an automatic differentiation tool developed by researchers at
[INRIA](https://team.inria.fr/ecuador/en/tapenade/) (the French National
Institute for Research in Computer Science and Automation).

Tapenade is designed to automatically generate derivative code for numerical
simulation programs written in Fortran or C. It enables the computation of
gradients, Hessians, and other derivatives efficiently, which is particularly
useful in fields such as optimization, sensitivity analysis, and scientific
computing.

By analyzing the source code of the original program, Tapenade generates code
that computes the derivatives of the program's outputs with respect to its
inputs. This capability is crucial in many scientific and engineering
applications where the ability to efficiently compute derivatives is essential.

Overall, Tapenade simplifies the process of incorporating automatic
differentiation into existing numerical simulation codes, making it a valuable
tool for researchers and developers working in computational science and
engineering.

## How we use it
In `yaeos` we developed a wrapper object that receives a set of routines from
a differentiated module and uses and internal logic to get the desired $A_r$ 
derivatives.

## Obtain a tapenade differentiated EoS
Getting a usable $A_r$ equation of state with `tapenade` is fairly easy.

1. Install `tapenade` and assure that you have the `tapenade` executable in
   your `PATH`.
2. Setup your new model following the [template file](./template.f90).
   a full implementation of the PengRobinson EoS can be seen at
   [pr.f90](./pr.f90) as an example.
3. Run the script `gen_tapemodel.sh`, providing your file as an argument:
   ```bash
   bash gen_tapemodel.sh <your_model_file.f90>
   ```
   This will generate a new folder `tapeout`, with your differentiated model
   inside.
4. Some little post-process must be done due to some details in the `tapenade`
   implementation. These are described in the base template but can also
   be checked on the [differentiated PR76 result after fixing the last details](./tapeout/pr_diff.f90)

To add your new tapenade model just include the file in your `src` folder and
use it with

```fortran
use yaeos, only: ArModel, pressure
use your_module_name, only: setup_model

class(ArModel), allocatable :: model

model = setup_model(<your parameters>)

call pressure(model, n, v, t)
```
