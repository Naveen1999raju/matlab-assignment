data = getgenbank('NC_045512');
features = featureparse(data);
featureview(data,{'gene'},'fontsize',12,'showpositions',true)
title('Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome.')
