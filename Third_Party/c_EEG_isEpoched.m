function isEpoched = c_EEG_isEpoched(EEG)
	if false
		isEpoched = EEG.trials > 1;
	else
		% handle special case where data was epoched to a single epoch (e.g. by trial rejection)
		% (infer that if there is an event exactly at time 0 and there is data available for negative times,
		%  then this was probably a single epoch)
		isEpoched = EEG.trials > 1 || (EEG.xmin < 0 && any(EEG.times(round([EEG.event.latency]))==0)); 
	end
end