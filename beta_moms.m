function Emn = beta_moms(m,n,a_ba,b_ba,b_cb)

Emn = (beta(m+n+a_ba,b_ba).*beta(n+a_ba+b_ba,b_cb))./...
(beta(a_ba,b_ba).*beta(a_ba+b_ba,b_cb));


end