function varargout = PropsSI(varargin)
  [varargout{1:max(1,nargout)}] = CoolProp_wrap(316,varargin{:});
end
