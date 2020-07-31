function varargout = get_input_pair_short_desc(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(156,varargin{:});
end
