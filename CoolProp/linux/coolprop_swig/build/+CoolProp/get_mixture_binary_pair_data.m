function varargout = get_mixture_binary_pair_data(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(159,varargin{:});
end
