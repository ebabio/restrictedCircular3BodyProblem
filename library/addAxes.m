function addAxes( axisScale )
% addAxes to any plot

plot3([0, axisScale, nan, 0, 0, nan, 0, 0], [0, 0, nan, 0, axisScale, nan, 0, 0], [0, 0, nan, 0, 0, nan, 0, axisScale], 'Color',[.2 .2 .2] );
text([axisScale, 0, 0], [0, axisScale, 0], [0, 0, axisScale], ['X';'Y';'Z']);

end

