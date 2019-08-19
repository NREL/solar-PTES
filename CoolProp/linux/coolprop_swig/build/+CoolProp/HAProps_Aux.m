function varargout = HAProps_Aux(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(425,varargin{:});
end
