function diff,arr
  n=n_elements(arr)
  delta=[0,arr(1:n-1)-arr(0:n-2)]
  delta[0]=delta[1]
  return,delta
end

;example of use if x is an array,
;delx=diff(x) gives differences between values
