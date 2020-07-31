function varargout = cair_sat(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(432,varargin{:});
end
