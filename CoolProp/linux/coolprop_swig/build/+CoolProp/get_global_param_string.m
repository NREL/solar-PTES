function varargout = get_global_param_string(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(323,varargin{:});
end
