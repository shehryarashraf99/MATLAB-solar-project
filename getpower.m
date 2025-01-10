function power = getpower(input,ambient_temp,panels_per_tracker,NOCT,Power_Coeff,STC_cell_temp,STC_eff,STC_OC,panel_area)
% input is the input irradiance in watts per meter square
% ambient temp is the site temperature
% panels per tracker is the number of panels in  the tracker half
if input==0
    power=0;
else


n=input/1000; % input irradiance
% STC_eff=0.192; % from panel specifications, panel efficiency at standard irradiance and temperature
% STC_temp=25; % standrd cell temperature
% Voc=44.4;% open circuit voltage at standard temperature
% NOCT=45; % nominal operating cell temperature

area= panel_area*10^-6;% this is the area of a panel in the tracker, in meter square .



%power_coefficient=-0.36; % from panel specifications

module_temp=ambient_temp+ (NOCT-STC_cell_temp)*n; % compute actual module temperature from irradiance and ambient temperature


B=1.38064852*10^-23; % boltzman's coeficient
q=1.60217662*10^-19; % charge of an electron
eff=STC_eff*(1- ((273+module_temp)*B*log(n)/(q*STC_OC))); % actual efficiency

Power_output=input*area*eff; % power output under current light intensity
loss=Power_Coeff/100 *(module_temp-STC_cell_temp)*Power_output; % power loss due to temperature
Power_output=Power_output+loss; % adjusted power output.


power=Power_output*panels_per_tracker; % total power output.
end
end