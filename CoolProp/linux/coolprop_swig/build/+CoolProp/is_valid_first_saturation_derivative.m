function varargout = is_valid_first_saturation_derivative(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(151,varargin{:});
end
