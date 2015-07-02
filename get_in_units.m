function val = get_in_units(objH,propName,unit)
%GET_IN_UNITS Get a graphics object property in the given unit
%   Detailed explanation goes here

oldUnits = get(objH,'Units'); 
set(objH,'Units',unit);   
val = get(objH,propName);   
set(objH,'Units',oldUnits);    

end

