printf("lib_cubic_surface requires:\n"):
printf("poly_helps.mpl\n\n"):

#assuming lines at infinity or those with x=0 does not occur!
find_gb_lines := proc(f,var::symbol:=x)
  local Lines, ideal:
  global x,y,z,a,b,c,d:
  if var = x then
    Lines := subs([y=a+b*x,z=c+d*x],f):
  elif var = z then #use to get constant x solutions
    Lines := subs([x=c+d*z,y=a+b*z],f):
  fi:
  ideal := [seq(coeff(Lines,var,i),i=0..3)]:
  Groebner:-Basis(ideal, plex(a,b,c,d));  
end proc:

#transforms a cubic form into one in general position (not coaxial with x,y,z and no lines at infinity)
generic_cubic_transform := proc(f,M, homog::boolean:=false)
  local g, G, newV:
  global w,x,y,z:
  
  if homog then g := Groebner:-Homogenize(f,w):
  else g:=f: fi:
  
  newV := M.<w,x,y,z>:
  g := subs([w = newV[1], x=newV[2], y=newV[3], z=newV[4]], g):
  G:=g:
  g := subs(w=1, g):
  g, G, find_gb_lines(g);
end proc:
