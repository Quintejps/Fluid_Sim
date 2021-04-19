# Fluid_Sim

Tasks

1. Add the capability to the simulator to create multiple sources for velocity and density without
modifying the source code. This means that you can either receive them via command line or via a
file. Either is accepted, you DON’T have to do both.

For this task, i decided to implement a file reader, which i prefer because is more simple for the user to give a bigger input to the fluid simulator. Python has a function called open() which helps you open files and save them into variables, then i used a for loop to read the input file line by line to get all the data for the new velocity, density or figure.

2. Animate the velocity forces. This should also be configurable via a command line parameter or can
be included in the input file. I would expect at least 2 behaviors.

For this task, i used math library to get the sin and cos, this helped me modify velocity equation.

3. Create color schemas for the simulation. Right now it is in yellow and blue, but I would want you
to have multiple coloring options.

For this task, i created a new cmap using RGB and then with help of a function registered the cmap to apply it into the imshow function. 

4. Simulate the presence of objects in the simulation. As before, this could be provided via command
line parameter or via file.

For this task, i created a for loop which loops over the figure area established in the input file and then setting the velocity and density functions equal to zero in each x and y.

Instructions for the input.txt file

1V-Fisrt velocity animation
1V(x,y)=(here goes x,here goes y)|velocity=(here goes velocity x,here goes velocity x)
Example:
1V(x,y)=(10,45)|velocity=(5,2)

2V-Second velocity animation
2V(x,y)=(here goes x,here goes y)|velocity=(here goes velocity x,here goes velocity x)
Example:
2V(x,y)=(10,45)|velocity=(5,2)

3V-Third velocity animation
3V(x,y)=(here goes x,here goes y)|velocity=(here goes velocity x,here goes velocity x)
Example:
3V(x,y)=(10,45)|velocity=(5,2)

D-Density
D(x1:x2,y1:y2)=(here goes x1:here goes x2,here goes y1:here goes y2)|density=here goes the density
Example:
D(x1:x2,y1:y2)=(0:20,40:60)|density=100

F-Figure
F(x1:x2,y1:y2)=(here goes x1:here goes x2,here goes y1:here goes y2)
Example:
F(x1:x2,y1:y2)=(15:45,15:45)



