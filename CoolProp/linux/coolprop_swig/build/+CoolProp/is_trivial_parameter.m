function varargout = is_trivial_parameter(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(148,varargin{:});
end
