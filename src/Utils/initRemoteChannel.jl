export initRemoteChannel

"""
function ref = initRemoteChannel(f,pid,args...)

Runs the function f with arguments args on worker

input: 
       func::Union{Function,Type} -- function or type constructor
        pid::Int64 -- worker that will run f and store result
       args::varargs -- Arguments of f. These must be available on worker pid
       kwargs::varargs -- Key word arguments of f. These must be available on worker pid
       
output:
        ref::RemoteChannel -- Reusable reference to output of f.
       
"""
function initRemoteChannel(func::Union{Function,Type}, pid::Int64, args...; kwargs...)
  return RemoteChannel(()->initChannel(func,args,kwargs), pid)
end

function initChannel(func::Union{Function,Type},args::Tuple,kwargs::Array)
  obj = func(args...; kwargs...)
  chan = Channel{typeof(obj)}(1)
  put!(chan,obj)
  return chan
end
