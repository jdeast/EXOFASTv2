function tcb2tdb, tcb

lb = 1.550519768d-8
tdb0 = -6.55d-5
t0 = 2443144.5003725d0

tdb = tcb - (lb*(tcb-t0)*86400 + tdb0)/86400

return, tdb

end
