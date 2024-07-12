module GeometricControl

using LinearAlgebra

export geometric_control, actuator_allocation

function hat_map(v::Vector)
    return [0.0 -v[3] v[2];
            v[3] 0.0 -v[1];
            -v[2] v[1] 0.0]
end
function vee_map(skew::Matrix)
    vec = 1/2 .* [skew[3,2]-skew[2,3];skew[1,3]-skew[3,1];skew[2,1]-skew[1,2]]
    return vec
end
function dir_cosine(q)
    C_B_I = zeros(3,3)
    C_B_I[1,1] = 1 - 2 * (q[3]^2 + q[4]^2)
    C_B_I[1,2] = 2 * (q[2] * q[3] + q[1] * q[4])
    C_B_I[1,3] = 2 * (q[2] * q[4] - q[1] * q[3])
    C_B_I[2,1] = 2 * (q[2] * q[3] - q[1] * q[4]) 
    C_B_I[2,2] = 1 - 2 * (q[2]^2 + q[4]^2) 
    C_B_I[2,3] = 2 * (q[3] * q[4] + q[1] * q[2])
    C_B_I[3,1] = 2 * (q[2] * q[4] + q[1] * q[3])
    C_B_I[3,2] = 2 * (q[3] * q[4] - q[1] * q[2])
    C_B_I[3,3] =  1 - 2 * (q[2]^2 + q[3]^2)
    return C_B_I'
end
function controller(t::Float64,x::Vector,kR::Float64,komega::Float64,kx::Float64,kv::Float64,param::Dict)
    g = param["g"]
    m = param["m"]
    J = param["J"]

    rx = x[1]
    ry = x[2]
    rz = x[3]

    quat = x[4:7]

    vx = x[8]
    vy = x[9]
    vz = x[10]

    p = x[11]
    q = x[12]
    r = x[13]
    omega = [p;q;r]

    # R
    # R = get_R([phi,theta,psi])
    R = dir_cosine(quat)

    # desired position (hard coding)
    frequency = 2*pi/5
    radius = 5

    b = 1 / 10
    zdes = 5.0
    zdes_dot = 0.0
    zdes_ddot = 0.0

    xdes = radius * cos(frequency*t)
    xdes_dot = - radius * frequency * sin(frequency*t)
    xdes_ddot = - radius * frequency^2 * cos(frequency*t)
    ydes = radius * sin(frequency*t)
    ydes_dot = radius * frequency * cos(frequency*t)
    ydes_ddot = - radius * frequency^2 * sin(frequency*t)

    psi_des = deg2rad(0)

    # errors
    ex = [rx;ry;rz] - [xdes;ydes;zdes]
    ev = [vx;vy;vz] - [xdes_dot;ydes_dot;zdes_dot]



    # total thrust
    xddot_des = [xdes_ddot;ydes_ddot;zdes_ddot]
    Fd = - kx * ex - kv * ev + m.*g.*[0;0;1] + m * xddot_des
    Fz = dot(Fd,R * [0;0;1])

    # desired attitudes
    if norm(Fd) < 1e-8
        @error("Fd is too small")
    end
    b3d = Fd ./ norm(Fd)
    b1d = [cos(psi_des);sin(psi_des);0]
    b2d = cross(b3d,b1d)
    b2d = b2d ./ norm(b2d)
    Rd = zeros(3,3)
    Rd[:,1] = cross(b2d,b3d)
    Rd[:,2] = b2d
    Rd[:,3] = b3d
    Rd_dot = dir_cosine([1;0;0;0])
    omega_des = zeros(3)
    omega_dot_des = zeros(3)

    # errors
    eR = 0.5 * vee_map(Rd'*R - R'*Rd)
    e_omega = omega - R'*Rd*omega_des

    M = (-kR .* eR - komega .* e_omega) + cross(omega,J*omega) - J*(hat_map(omega)*R'*Rd*omega_des - R'*Rd*omega_dot_des)

    return Fz,M
end


function actuator_allocation(Fz::Float64,M::Vector,N_inv::Matrix)
    return N_inv * [Fz;M]
end

end