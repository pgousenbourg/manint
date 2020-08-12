Model3D

Author: Steven Michael (smichael@ll.mit.edu)
Date:   5/19/05

Description

Copyright 2005 MIT Lincoln Laboratory 
Released under the GNU LGPL.  See lgpl.txt in this directory


The following MATLAB class provides the ability to load 3D studio MAX (.3ds)
and AutoCad (DXF) files into MATLAB.  Routines exist for plotting the
resulting class using OpenGL and for very simple data manipulation (rotation,
translation, superposition of models, etc...)


Usage:

The model can be loaded using the 'model3d' command.  Allowed manipulations of
the model are described in the accompanying '.m' files.  

Example:

%Load a 3ds model 'sample.3ds' and plot it,
m = model3d('sample.3ds')

m =

        model3d object: 1-by-1

>> plot(m3)
% Create a new model shifted by [10 10 10] in X,Y,Z
>> m2 = m+[10 10 10];
% Superimpose the shifted model on the original
% and output it in m2
>> m2 = m+m2;
% Magnify the model by 3
>> m2 = magnify(m,3);
% Also magnifies the model by 3
>> m2 = m * 3;
% Rotate the model by PI about the Z axis
>> m2 = qrot(m,[0 0 1],pi);
    
Note: 

The Windows files with the ".mexw32" extension may have to be 
renamed to a ".dll" extension to work with versions of MATLAB older than
R14 service pack 3.

Note: 

Under windows, the files are compiled with the Microsoft Visual 
Studio 8.0 compilers.  This may require installation of the Visual Studio
runtime libraries. These are included with MATLAB under:
$MATLAB/bin/win32/vcredist_x86.exe
See: http://www.mathworks.com/support/solutions/data/1-2223MW.html 

