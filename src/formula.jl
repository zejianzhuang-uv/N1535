

"""redefine sqrt so that its cut is along the positive x axis"""
function xsqrt(x)
    imag(x) >=0 ? sqrt(x+0im) : -sqrt(x-0im)#imag(x) >=0 ? sqrt(x+0im) : -sqrt(x-0im)
end


"""Kallen function"""
function kallen(x, y, z)
    return x^2 + y^2 + z^2 - 2x*y - 2y*z - 2z*x
end

"""energy of baryon"""
function baryon_energy(w, mchi)
    m, M = mchi
    return (w^2 + M^2 - m^2) / (2w)
end

"""energy of meson"""
function meson_energy(w, mchi)
    m, M = mchi
    return (w^2 + m^2 - M^2) / (2w)
end


"""memuntum in center mass"""
function qcm(w, m1, m2)
    return xsqrt(kallen(w^2, m1^2, m2^2) ) / (2w)
end




""" Convert a complex number to a tuple """
function cmplx(xx)
    return real(xx), imag(xx)
end