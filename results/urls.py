from django.contrib import admin
from django.urls import include, re_path
from django.conf import settings
from django.conf.urls.static import static
from .views import Results,queryPlotHTML,queryPlotFeature,queryPlotTopDE,Results_tutorial

urlpatterns = [
    re_path(r'^tutorial2$', Results_tutorial.as_view(), name='tutorial_2'),
    re_path(r'^ajax_plot_topDE$', queryPlotTopDE, name='ajax_plot_topDE'),
    re_path(r'^ajax_plot$', queryPlotHTML, name='ajax_plot'),
    re_path(r'^ajax_feature$', queryPlotFeature, name='ajax_feature'),
    re_path('', Results.as_view(), name='results'),
    ]+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)