module ContingencyTableCorectnessTest

using DataFrames
using ContingencyTables
using Base.Test

  
  #
  #	Concrete Testing Data 
  #
  #
  
  A = ContingencyTable(DataFrame(quote
  a = [0,1]
  probab = [0.2,0.8]
  end))
  
  A2 = ContingencyTable(DataFrame(quote
  a = ["A","B"]
  probab = [0.2,0.8]
  end))
  
  AB = ContingencyTable(DataFrame(quote
  a = [0,0,1,1]
  b = [0,1,0,1]
  probab = [0.06,0.14,0.24,0.56]
  end))
  
  A_E = ContingencyTable(DataFrame(quote
  a = [0,1]
  probab = [0.2,0.8]
  end),[("b",x->x==0)])
  
  A_BE = ContingencyTable(DataFrame(quote
  a = [0,1,0,1]
  b = [0,0,1,1]
  probab = [0.1,0.4,0.2,0.3]
  end))
  
  A_EC = ContingencyTable(DataFrame(quote
  a = [0,0,0,0,1,1,1,1]
  c = [0,0,1,1,0,0,1,1]
  d = [0,1,0,1,0,1,0,1]
  probab = [0.04,0.06,0.08,0.02,0.16,0.24,0.32,0.08]
  end),[("b",x->x==0)])

  A_E2 = ContingencyTable(DataFrame(quote
  a = [0,1]
  probab = [0.2,0.8]
  end),[("c",x->x==0)])
  
  A_E3 = ContingencyTable(DataFrame(quote
  a = [0,1]
  probab = [0.2,0.8]
  end),[("c",x->x==2)])
  
  A_E4 = ContingencyTable(DataFrame(quote
  a = [0,1]
  probab = [0.2,0.8]
  end),[("c",x->x+2==2)])
  
  A_E5 = ContingencyTable(DataFrame(quote
  a = [0,1]
  probab = [0.2,0.8]
  end),[("b",x->x==0),("c",x->x==1)])
  
  A_E5_clone = ContingencyTable(DataFrame(quote
  a = [0,1]
  probab = [0.2,0.8]
  end),[("b",x->x==0),("c",x->x==1)])
  
  A_E6 = ContingencyTable(DataFrame(quote
  a = [0,1]
  probab = [0.2,0.8]
  end),[("c",x->x==1),("b",x->x==0)])

  A_E7 = ContingencyTable(DataFrame(quote
  a = [0,1]
  probab = [0.6,0.4]
  end),[("c",x->x==1),("b",x->x==0)])
  
  A_BE7 = ContingencyTable(DataFrame(quote
  a = [0,1,0,1]
  b = [0,0,1,1]
  probab = [0.15/0.375,0.1/0.375,0.05/0.375,0.075/0.375]
  end),[("c",x->x==1)])
  
  A_CE7 = ContingencyTable(DataFrame(quote
  a = [0,1,0,1]
  c = [0,0,1,1]
  probab = [0.2/0.7,0.25/0.7,0.15/0.7,0.1/0.7]
  end),[("b",x->x==0)])

  A_BCE7 = ContingencyTable(DataFrame(quote
  a = [0,0,0,0,1,1,1,1]
  b = [0,0,1,1,0,0,1,1]
  c = [0,1,0,1,0,1,0,1]
  probab = [0.2,0.15,0.1,0.05,0.25,0.1,0.075,0.075]
  end))  
  
  B = ContingencyTable(DataFrame(quote
  b = [0,1]
  probab = [0.3,0.7]
  end))
  
  BC_C = ContingencyTable(DataFrame(quote
  b = [0,0,0,0,1,1,1,1]
  c = [0,0,1,1,0,0,1,1]
  d = [0,1,0,1,0,1,0,1]
  probab = [0.1,0.225,0.2,0.075,7/30,0.525,7/15,0.175]
  end),Array((String,Function),0),["d"])
  
  C = ContingencyTable(DataFrame(quote
  c = [0,0,1,1]
  d = [0,1,0,1]
  probab = [0.2,0.3,0.4,0.1]
  end))
  
  C_C = ContingencyTable(DataFrame(quote
  c = [0,0,1,1]
  d = [0,1,0,1]
  probab = [1/3,0.75,2/3,0.25]
  end),Array((String,Function),0),["d"])
  
  D = ContingencyTable(DataFrame(quote
  d = [0,1,0,1]
  c = [0,0,1,1]
  probab = [0.2,0.3,0.4,0.1]
  end))
  
  E = ContingencyTable(DataFrame(quote
  c = [1,1,0,0]
  d = [0,1,0,1]
  probab = [0.4,0.1,0.2,0.3]
  end))
  
  E_C = ContingencyTable(DataFrame(quote
  c = [1,1,0,0]
  d = [0,1,0,1]
  probab = [0.4,0.1,0.2,0.3]
  end),Array((String,Function),0),["d"])
  
  E_C_clone = ContingencyTable(DataFrame(quote
  c = [1,1,0,0]
  d = [0,1,0,1]
  probab = [0.4,0.1,0.2,0.3]
  end),Array((String,Function),0),["d"])
  
  E_C2 = ContingencyTable(DataFrame(quote
  c = [1,1,0,0]
  d = [0,1,0,1]
  probab = [0.4,0.1,0.2,0.3]
  end),Array((String,Function),0),["c"])
  
  E_C3 = ContingencyTable(DataFrame(quote
  c = [1,1,0,0]
  d = [0,1,0,1]
  probab = [0.4,0.1,0.2,0.3]
  end),Array((String,Function),0),["c","d"])
  
  E_C4 = ContingencyTable(DataFrame(quote
  c = [1,1,0,0]
  d = [0,1,0,1]
  probab = [0.4,0.1,0.2,0.3]
  end),Array((String,Function),0),["d","c"])
    
  handler(r::Success) = println("Success on $(r.expr)")
  handler(r::Failure) = println(r.expr),error("Error on custom handler: $(r.expr)")
  handler(r::Error)   = rethrow(r)
  
  withhandler(handler) do
    #
    #	testing isequal function
    #
	@test isequal(A,A)	#ConTables are themselves
	@test isequal(B,B)
	@test !isequal(A,B)	#different tables are different
	@test !isequal(B,A)
	@test !isequal(A,A2)
	@test !isequal(A,A_E)	#with/without evidence
	@test !isequal(E,E_C)	#unconditioned/conditioned
	@test isequal(C,D)	#tables with different row order are equal
	@test isequal(D,C)
	@test isequal(C,E)	#tables with different line order are equal
	@test isequal(E,C)
	@test isequal(D,E)	#tables with different row and line order
	@test isequal(E,D)	#are equal
	@test isequal(A_E,A_E)	#ConTables with evidence are themselves
	@test isequal(A_E5,A_E5)
	@test isequal(A_E5,A_E5_clone)	#this first caused problems, as in Functions.code.ast the line the function was declared was stored as well 
					#and therefore my comparison would consider them as different, fixed by excluding that value in comparison
	@test isequal(A_E5,A_E6)#with different evidence order
	@test !isequal(A_E,A_E2)#ConTables with different evidence are different
	@test !isequal(A_E,A_E3)#different evidence variable
	@test !isequal(A_E2,A_E3)#different evidence function
	@test !isequal(A_E2,A_E4) #note: functions with same output but different syntax are considered different!
	@test !isequal(A_E,A_E5)#more evidence vs less evidence
	@test isequal(E_C,E_C)	#conditioned ConTables are themselves
	@test isequal(E_C_clone,E_C)#this was neccessary to see whether the problem of lines 145-146 exists for conditions as well
	@test isequal(E_C3,E_C3)
	@test !isequal(E_C,E_C2)#tables with different conditions are different
	@test !isequal(E_C2,E_C3)#condition is subset of 2nd tables conditions
	@test isequal(E_C3,E_C4)#tables with different condition order are same
	
    #
    #	in the following tests
    #	fuzzy is set to true, if neccessary,
    #	to allow small variance in probabilities
    #	(equal to constant EPSILON in ContingencyTable.jl)
    #	necessary here means the test was unsuccessful otherwise
    #
    #
    #	testing fproduct
    #
	@test isequal(fproduct(A,B),AB,true)		#correctness
	@test isequal(fproduct(B,C_C),BC_C,true)	#including conditions
	@test isequal(fproduct(A_E,C),A_EC,true)	#including evidence
	@test isequal(fproduct(A,B),fproduct(B,A))	#symmetry
	@test isequal(fproduct(A,E_C),fproduct(E_C,A))	#symmetry including conditions
	@test isequal(fproduct(C,A_E5),fproduct(A_E5,C))	#symmetry including evidence
	@test isequal(fproduct(fproduct(A,B),C),fproduct(fproduct(B,C),A),true)#transitivity
	#@test isequal(fproduct(fproduct(A,B),A),fproduct(A,B))	#absorption now only works for conditional tables
	@test !isequal(fproduct(fproduct(A,C),fproduct(A,B)),fproduct(fproduct(A,B),C),true)#overlapping does make a difference
    
    #
    #	testing condition
    #
    #
	@test isequal(C_C,condition(C,["d"]),true)	#correctness
	@test isequal(condition(condition(fproduct(A,C),["c"]),["d"]),condition(condition(fproduct(A,C),["d"]),["c"]),true)#simmetry
	@test isequal(condition(fproduct(A,C),["c","d"]),condition(condition(fproduct(A,C),["c"]),["d"]),true) #at once vs iterative
	
    #
    #	testing evidence
    #
    #
	@test isequal(A_E,evidence(A_BE, [("b",x->x==0)]))	#correctness
	@test isequal(evidence(A_BCE7, [("b",x->x==0),("c",x->x==1)]),A_E7,true)
	@test isequal(evidence(evidence(A_BCE7,[("b",x->x==0)]),[("c",x->x==1)]),evidence(evidence(A_BCE7,[("c",x->x==1)]),[("b",x->x==0)]),true)#simmetry
	@test isequal(evidence(A_BCE7, [("b",x->x==0),("c",x->x==1)]),evidence(evidence(A_BCE7,[("b",x->x==0)]),[("c",x->x==1)]),true)#at once vs iterative
    
    #
    #	testing marginalize
    #
    #
	@test isequal(A, marginalize(fproduct(A,B),"b")) #correctness
	@test isequal(B, marginalize(fproduct(A,B),"a"),true)
	@test isequal(AB, fproduct(marginalize(AB,"a"),marginalize(AB,"b")))#note: this only works for independent random variables
	@test !isequal(C, fproduct(marginalize(C,"c"),marginalize(C,"d")),true)#see above
	@test isequal(marginalize(marginalize(fproduct(AB,C),"a"),"c"),marginalize(marginalize(fproduct(AB,C),"c"),"a"),true)#symmetry
	
    #
    #	testing some correlations
    #
    #
	@test isequal(marginalize(condition(fproduct(A,B),["b"]),"b"),A,true)
  end
end