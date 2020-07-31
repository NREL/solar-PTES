function varargout = HAPropsSI(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(423,varargin{:});
end
