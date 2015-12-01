function output_txt = NormBarCallback(obj,event_obj,slbls,sidx)

pos = get(event_obj,'Position');
output_txt=slbls{sidx(pos(1))};
get(event_obj);