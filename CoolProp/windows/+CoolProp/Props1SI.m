function varargout = Props1SI(varargin)
  [varargout{1:max(1,nargout)}] = CoolPropMATLAB_wrap(344,varargin{:});
end
