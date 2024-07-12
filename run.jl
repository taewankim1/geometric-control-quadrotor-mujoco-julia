using Plots
# using LinearAlgebra
using MuJoCo
import MuJoCo: step! as mj_step

include("./src/GeometricControl.jl")

model = load_model("./quadrotor.xml")
data = init_data(model)

forward!(model,data)

# Parameters
H = mj_zeros(model.nv, model.nv)
L = mj_fullM(model,H , data.qM)
param = Dict()
param["J"] = H[end-2:end,end-2:end]
param["m"] = 1.325
param["g"] = 9.81

# Parameters for actuators
dx = 0.14
dy = 0.18
c = 0.0201
N = [1 1 1 1;-dy dy dy -dy;dx dx -dx -dx;-c c -c c]
N_inv = inv(N)

# Control gains
kR = 30.0 
komega = 5.0

kx = 5.0
kv = 2.0

dt_control = model.opt.timestep
sim_time = 30.0
t_ctrl = range(0,sim_time,step=dt_control)

reset!(model,data)
ctrl_states = zeros(13, length(t_ctrl))
for (idx,t_) in enumerate(t_ctrl)
    @assert(isapprox(t_, data.time))
    ctrl_states[:,idx] = get_physics_state(model, data)
    Fz, M = GeometricControl.controller(t_,ctrl_states[:,idx],kR,komega,kx,kv,param)
    data.ctrl .= GeometricControl.actuator_allocation(Fz,M,N_inv)
    mj_step(model, data)
end

plot(ctrl_states[1,:],ctrl_states[2,:],title="position")
xlabel!("X (m)")
ylabel!("Y (m)")

init_visualiser()
visualise!(model, data, trajectories = ctrl_states)