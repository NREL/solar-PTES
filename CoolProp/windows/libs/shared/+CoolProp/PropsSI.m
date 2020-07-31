function varargout = PropsSI(varargin)
  [varargout{1:max(1,nargout)}] = CoolPropMATLAB_wrap(345,varargin{:});
end
