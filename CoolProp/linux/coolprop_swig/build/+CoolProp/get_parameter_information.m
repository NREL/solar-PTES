function varargout = get_parameter_information(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(145,varargin{:});
end