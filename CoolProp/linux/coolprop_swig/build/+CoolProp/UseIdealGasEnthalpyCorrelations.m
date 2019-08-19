function varargout = UseIdealGasEnthalpyCorrelations(varargin)
  [varargout{1:nargout}] = CoolProp_wrap(429,varargin{:});
end
