

include("formula.jl")
include("LECs_charged_basis.jl")

""" Potential up to next to leading order (S-wave) """
function V_up_to_nlo_swave(w, par, C, D, L; n=2)
    wt = WT(w, par, C, n=n)
    vnlo = vnlo_swave(w, par, D, L, n=n)
    vb = Born_s(w, par, n=n) + Born_u(w, par, n=n)
    return wt + vnlo + vb
end

function V_lo(w, par, lec; n=2)
    vb = Born_s(w, par, n=n) + Born_u(w, par, n=n)
    wt = WT(w, par, lec, n=n)
    return vb + wt
end

function WT(w, par, lec; n=2)
    
    #lec = lec_lo()
    f = par[:decons]
    mch = par[:mch]
    v = zeros(ComplexF64, n, n)
    for i in 1:n
        for j in 1:n
            Ni = norm_factor(w, mch[i] )
            Nj = norm_factor(w, mch[j] )
            v[i, j] = (2w - mch[i][2] - mch[j][2] ) * Ni * Nj * lec[i, j]
        end
    end
    
    return -1 / (4f^2) * v 
end

function vnlo_swave(w, par, D, L; n=2)
    f = par[:decons]
    mch = par[:mch]
    v = zeros(ComplexF64, n, n)
    for i in 1:n
        for j in 1:n
            Ni = norm_factor(w, mch[i] )
            Nj = norm_factor(w, mch[j] )
            #NNi = mch[i][2] + baryon_energy(w, mch[i])
            #NNj = mch[j][2] + baryon_energy(w, mch[j])

            ωi = meson_energy(w, mch[i])
            ωj = meson_energy(w, mch[j] )
            v[i, j] = Ni * Nj * (D[i, j] - 2*L[i, j]*ωi*ωj)
            #Ni * Nj * (D[i, j] - 2*L[i, j]*(ωi*ωj + qcm(w, mch[i]...)^2 * qcm(w, mch[j]...)^2 / (3NNi*NNj) ) )
            # The contributation from the term qi^2 qj^2 / (2Na*Nb) can be ignored.
        end
    end
    return v * 1/f^2
end

function Born_s(w::Any, params::Dict; n=2)
    f = params[:decons]
    mch = params[:mch]
    mB = params[:mB]
    v = zeros(ComplexF64, n, n)
    for i in 1:n
        for j in 1:n
            Ni = norm_factor(w, mch[i])
            Nj = norm_factor(w, mch[j])
            Mi = mch[i][2]
            Mj = mch[j][2]
            v[i, j] = Ni *Nj *(w - Mi) *(w - Mj) * sum([LECs_cBorn_s(i, j)[k] / (w + mB[k]) for k in 1:6])
        end
    end
    return v * (1. /(12e0f^2) )
end

function FFF1(s, mi, mj, Ei, Ej, qi, qj, mB)
    EPS = 8.
    if abs((s + mB^2 - mi^2 -mj^2 -2e0Ei*Ej + 2e0*qi*qj)) <= 1e-8
        s = (sqrt(s) - EPS)^2
        return log((s + mB^2 - mi^2 -mj^2 -2e0Ei*Ej - 2e0qi*qj) / (s + mB^2 - mi^2 -mj^2 -2e0Ei*Ej + 2e0*qi*qj) )
    else
        return log((s + mB^2 - mi^2 -mj^2 -2e0Ei*Ej - 2e0qi*qj) / (s + mB^2 - mi^2 -mj^2 -2e0Ei*Ej + 2e0*qi*qj) )
    end
    
end

function Born_u(w::Any, params::Dict; n=2)
    # println("The U channel")
    LECs_Born = LECs_cBorn_u
    EPS = 1e-7
    f = params[:decons]
    mch = params[:mch]
    mB = params[:mB]
    s = w*w
    v = zeros(ComplexF64, n, n)
    for i in 1:n
        for j in 1:n
            Ni = norm_factor(w, mch[i] )
            Nj = norm_factor(w, mch[j] )
            Ei = baryon_energy(w, mch[i] )
            Ej = baryon_energy(w, mch[j] )
            mi, Mi = mch[i]
            mj, Mj = mch[j]
            qi = (qcm(w, mch[i]...) )
            qj = (qcm(w, mch[j]...) )

            v[i, j] = sum([(w + mB[k] - (Mi + mB[k])*(Mj + mB[k])*(w + Mi + Mj - mB[k]) /(2e0*(Mi + Ei)*(Mj + Ej)) + (Mi + mB[k])*(Mj + mB[k]) /(4e0*qi*qj+EPS) *((w - Mi - Mj + mB[k]) - (s + mB[k]^2 -mi^2 -mj^2 -2e0Ei*Ej) /(2e0*(Mi + Ei)*(Mj + Ej)) *(w + Mi + Mj - mB[k])) * 
            FFF1(s, mi, mj, Ei, Ej, qi, qj, mB[k]) ) * LECs_Born(i, j)[k] for k in 1:6]) * Ni * Nj
        
        end
    end
    return (-1. / (12e0f^2)) * v
end

function norm_factor(w, mchi)
    EB = baryon_energy(w, mchi)
    M = mchi[2]
    
    return sqrt((EB + M) / 2M)
end


function lecs_cij()
    lec = [
        [2, 1, 1 / 2, 1, 0, sqrt(3) / 2],
        [1, 2, 2, 0, 0, 0],
        [1 / 2, 2, 0, 1 / 2, 2, 0],
        [1, 0, 1 / 2, 2, 1, -sqrt(3) / 2],
        [0, 0, 2, 1, 2, 0],
        [sqrt(3) / 2, 0, 0, -sqrt(3) / 2, 0, 0]]
    return vcat(lec'...) # convert the vector to a 6 x 6 matrix
end

#=
function lec_bi(par, b::Tuple)
    (mπ, mK) = par[:meson_mass]
    (b0, bD, bF) = b
    mu1 = sqrt(mK^2 + mπ^2)
    #mu2 = sqrt(5mK^2 - 3mπ^2)
    bb = [
        [4(b0 + bD)*mK^2, (bD - bF)*mu1^2, (bD - bF)*mu1^2 / 2, 2(bD + bF)*mK^2, 0, -(bD + 3bF)*mu1^2 / (2sqrt(3) )],
        [(bD - bF) * mu1^2, 4(b0 + bD)*mπ^2, 0, 0, 0, 0],
        [(bD - bF)*mu1^2 / 2, 0, 4(b0+bD)mπ^2, (bD - bF)*mu1^2/2, 0, 0],
        [2(bD + bF)*mK^2, 0, (bD - bF)*mu1^2/2, 4(b0 + bD)*mK^2, (bD - bF)*mu1^2, (bD + 3bF)*mu1^2/(2sqrt(3))],
        [0, 0, 0, (bD - bF)*mu1^2, 4(b0+bD)*mπ^2, 0],
        [-(bD + 3bF)*mu1^2 / (2sqrt(3) ), 0, 0, (bD + 3bF)*mu1^2 / (2sqrt(3)), 0, 4(3b0+bD)mπ^2 / 3]
    ]
    return vcat(bb'...) #hcat(bb...)
end
=#
function lec_bi(par, b)
    (mπp, mπ0, mπm, mKbar0, mKm) = par[:meson_mass]
    (b0, bD, bF) = b .* 1e-3
    mu1(mK, mπ) = sqrt(mK^2 + mπ^2)
    #mu2 = sqrt(5mK^2 - 3mπ^2)
    bb = [
        [4(b0 + bD)*mKm^2, (bD - bF)*mu1(mKm, mπm)^2, (bD - bF)*mu1(mKm, mπ0)^2 / 2, 2(bD + bF)*mKm^2, 0, -(bD + 3bF)*mu1(mKm, mπ0)^2 / (2sqrt(3) )],

        [(bD - bF) * mu1(mKm, mπm)^2, 4(b0 + bD)*mπm^2, 0, 0, 0, 0],

        [(bD - bF)*mu1(mKm, mπ0)^2 / 2, 0, 4(b0+bD)mπ0^2, (bD - bF)*mu1(mKbar0, mπ0)^2/2, 0, 0],

        [2(bD + bF)*mKm^2, 0, (bD - bF)*mu1(mKbar0, mπ0)^2/2, 4(b0 + bD)*mKbar0^2, (bD - bF)*mu1(mKbar0, mπp)^2, (bD + 3bF)*mu1(mKbar0, mπ0)^2/(2sqrt(3))],

        [0, 0, 0, (bD - bF)*mu1(mKbar0, mπp)^2, 4(b0+bD)*mπp^2, 0],
        
        [-(bD + 3bF)*mu1(mKm, mπ0)^2 / (2sqrt(3) ), 0, 0, (bD + 3bF)*mu1(mKbar0, mπ0)^2 / (2sqrt(3)), 0, 4(3b0+bD)mπ0^2 / 3]
    ]
    return vcat(bb'...) #hcat(bb...)
end

function lec_di(d)
    (d1, d2, d3, d4) = d .* 1e-3
    lec = [
        [2d2+d3+2d4, -d1+d2+d3, (-d1-d2+2d3)/2, d1+d2+d3, -2d2+d3, -sqrt(3)*(d1+d2) / 2],

        [-d1+d2+d3, 2d2+d3+2d4, -2d2+d3, -2d2+d3, -4d2+2d3, 0],

        [(-d1-d2+2d3) / 2, -2d2 + d3, 2(d3+d4), (-d1-d2+2d3)/2, -2d2+d3, 0],

        [d1+d2+d3, -2d2+d3, (-d1-d2+2d3)/2, 2d2+d3+2d4, -d1+d2+d3, -sqrt(3)*(d1+d2)/2],

        [-2d2+d3, -4d2+2d3, -2d2+d3, -d1+d2+d3, 2d2+d3+2d4, 0],

        [-sqrt(3)*(d1+d2)/2, 0, 0, -sqrt(3)*(d1+d2)/2, 0, 2d4]
    ]
    return vcat(lec'...) #hcat(lec...)
end


function Dmat(par, b)
    (mπ, mK) = par[:meson_mass]
    (b0, bD, bF) = b .* 1e-3
    μ1 = sqrt(mK^2 + mπ^2)
    μ2 = sqrt(5mK^2 - 3mπ^2)
    μ3 = sqrt(4mK^2 - mπ^2)
    μ4 = sqrt(16mK^2 - 7mπ^2)

    D11 = 4 * (b0 + bD) * mπ^2
    D12 = -sqrt(3 / 2) * (bD - bF) * μ1^2
    D13 = -(4bD * mπ^2) / sqrt(3)
    D14 = sqrt(3 / 2) * (bD + bF) * μ1^2

    D21 = D12
    D22 = 2 * (2b0 + 3bD + bF) * mK^2
    D23 = (bD + 3bF) * μ2^2 / (3sqrt(2))
    D24 = 0e0

    D31, D32 = D13, D23
    D33 = (4 / 9) * (3b0 * μ3^2 + bD * μ4^2)
    D34 = -(bD - 3bF) * μ2^2 / (3sqrt(2))

    D41, D42, D43 = D14, D24, D34
    D44 = 2 * (2b0 + 3bD - bF) * mK^2

    D = [D11 D12 D13 D14; D21 D22 D23 D24; D31 D32 D33 D34; D41 D42 D43 D44]
    return D
end

function Lmat(d)
    d1, d2, d3, d4 = d .* 1e-3#d[:d1], d[:d2], d[:d3], d[:d4]

    L11 = -4d2 + 4d3 + 2d4
    L12 = sqrt(3 / 2) * (d1 + d2 - 2d3)
    L13 = -sqrt(3)d3
    L14 = sqrt(3 / 2) * (d1 - d2 + 2d3)

    L21 = L12

    L22 = d1 + 3d2 + 2 * (d3 + d4)
    L23 = (d1 - 3d2 + 2d3) / sqrt(2)
    L24 = (6d2 - 3d3)

    L31, L32 = L13, L23

    L33 = 2 * (d3 + d4)

    L34 = (d1 + 3d2 - 2d3) / sqrt(2)

    L41, L42, L43 = L14, L24, L34

    L44 = -d1 + 3d2 + 2 * (d3 + d4)

    L = [L11 L12 L13 L14; L21 L22 L23 L24; L31 L32 L33 L34; L41 L42 L43 L44]

    return L
end

function C_mat()
    c = [
        [4, -sqrt(3/2), 0, sqrt(3/2)],
        [-sqrt(3/2), 3, 3/sqrt(2), 0],
        [0, 3/sqrt(2), 0, -3/sqrt(2)],
        [sqrt(3/2), 0, -3/sqrt(2), 3]
    ]
    return hcat(c...)
end