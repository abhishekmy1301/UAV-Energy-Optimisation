%TRANSITIONUAVFIDELITYOPENAUTOPILOT opens the outer loop autopilot. This
%autopilot works for both the mid fidelity and high fidelity plant model.
modelname='Outer_Loop_Autopilot';
if ~bdIsLoaded(modelname)
    open_system(modelname);
end

